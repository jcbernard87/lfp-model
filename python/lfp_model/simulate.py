"""Time integration and output (docs/model.md sections 7 and 9).

`run()` reproduces the original constant-current discharge: one linearized solve
per step, output rows written before each step, and the original exit tests.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np

import bandsolver

from . import kinetics
from .model import Assembler, C, CS, P1, P2
from .params import Params, f32

HEADER = ("State", "Time", "Voltage", "Equivalence", "Anode_Eta", "anode_exchange_c", "Edge_c0")
UNITS = ("CDR", "hours", "Volts", "electron_equivs", "mV", "mA/cm2", "mol/cm3")


def _fmt_fixed(v: float) -> str:
    return f"{'NaN':>12}" if math.isnan(v) else f"{v:12.5f}"


def _fmt_sci(v: float) -> str:
    return f"{'NaN':>15}" if math.isnan(v) else f"{v:15.5E}"


def format_header() -> str:
    def row(cols):
        # Fortran A<w> right-justifies shorter strings and keeps the leftmost w characters of longer ones
        return (f"{cols[0][:5]:>5} " + " ".join(f"{c[:12]:>12}" for c in cols[1:3]) + " "
                + " ".join(f"{c[:15]:>15}" for c in cols[3:]))
    return row(HEADER) + "\n" + row(UNITS) + "\n"


def format_row(state: str, values) -> str:
    t, v, *rest = values
    return f"{state:>5} {_fmt_fixed(t)} {_fmt_fixed(v)} " + " ".join(_fmt_sci(x) for x in rest) + "\n"


@dataclass
class Result:
    rows: list = field(default_factory=list)      # (state, t_h, V, equiv, eta_mV, i0_mA, c_edge)
    exit_reason: str = ""
    steps: int = 0
    final_state: Optional[np.ndarray] = None

    @property
    def array(self) -> np.ndarray:
        return np.array([r[1:] for r in self.rows], dtype=float)

    def write(self, path) -> None:
        with open(path, "w") as fh:
            fh.write(format_header())
            for r in self.rows:
                fh.write(format_row(r[0], r[1:]))


def initial_state(p: Params) -> np.ndarray:
    c = np.empty((p.nj, 4))
    c[:, C] = p.c_bulk
    c[:, P1] = p.phi1_init
    c[:, P2] = p.phi2_init
    c[:, CS] = p.cs_init
    return c


def _output_row(p: Params, state: str, t: float, c: np.ndarray, mAhg: float):
    lit36 = f32(3.6) if p.mode == "faithful" else 3.6
    equiv = mAhg * p.M * lit36 / p.F
    i0_li = float(kinetics.li_exchange_current(p, c[0, C]))
    alpha = 0.5
    with np.errstate(invalid="ignore", divide="ignore"):
        if state == "C":
            eta = 0.5 * math.log(p.i_app / i0_li) / (alpha * p.F / (p.R * p.T))
        elif state == "D":
            eta = -(0.5 * math.log(p.i_app / i0_li)) / (alpha * p.F / (p.R * p.T))
        else:
            eta = 0.0
    return (state, t / float(3600), c[-1, P1] + eta, equiv, eta * 1.0e3, i0_li * 1.0e3, c[0, C])


class SolverFailure(RuntimeError):
    pass


def newton_step(asm: Assembler, c_old: np.ndarray, dt: float, *, backend: str = "fortran",
                pivot: str = "partial") -> tuple[np.ndarray, int]:
    """One backward-Euler step solved to convergence with Newton's method (corrected mode, fixes D-7)."""
    p = asm.p
    T = asm.time_terms(dt)
    scale = np.array([p.c_bulk, 1.0, 1.0, kinetics.cs_max(p)])
    c = c_old.copy()
    for k in range(1, p.newton_max_iter + 1):
        A, B, D, G, _ = asm.assemble(c, dt)
        G = G - T * (c - c_old)
        try:
            dc = bandsolver.solve(A, B, D, G, pivot=pivot, backend=backend)
        except (bandsolver.NonFiniteError, bandsolver.SingularBlockError) as e:
            raise SolverFailure(str(e)) from e
        lam = _bounded_step(c, dc, kinetics.cs_max(p))
        c = c + lam * dc
        if lam == 1.0 and np.max(np.abs(dc) / scale) <= p.newton_tol:
            return c, k
    raise SolverFailure(f"Newton did not converge in {p.newton_max_iter} iterations")


def _bounded_step(c, dc, csmax, keep=0.9):
    """Largest step length <= 1 keeping 0 < c, 0 < cs < cs_max (moves at most `keep` of the way to a bound)."""
    lam = 1.0
    for col, lo, hi in ((C, 0.0, None), (CS, 0.0, csmax)):
        x, d = c[:, col], dc[:, col]
        with np.errstate(divide="ignore", invalid="ignore"):
            neg = d < 0
            if np.any(neg):
                lam = min(lam, float(np.min(keep * (x[neg] - lo) / -d[neg])))
            if hi is not None:
                pos = d > 0
                if np.any(pos):
                    lam = min(lam, float(np.min(keep * (hi - x[pos]) / d[pos])))
    return max(lam, 0.0)


def advance(asm: Assembler, c: np.ndarray, dt: float, *, backend: str, pivot: str,
            stop=None, min_dt: float = 1.0e-6):
    """Advance by dt with backward Euler, halving the sub-step on Newton failure.

    `stop(c)` is checked after every sub-step; when it returns True the step ends early.
    Returns (c, time advanced, stopped).
    """
    t_done, h = 0.0, dt
    while t_done < dt:
        h = min(h, dt - t_done)
        try:
            c_new, _ = newton_step(asm, c, h, backend=backend, pivot=pivot)
        except SolverFailure:
            if h / 2 < min_dt:
                raise
            h = h / 2
            continue
        c, t_done = c_new, t_done + h
        if stop is not None and stop(c):
            return c, t_done, True
    return c, t_done, False


def run(p: Params, *, backend: str = "fortran", pivot: Optional[str] = None, max_steps: Optional[int] = None) -> Result:
    """Constant-current discharge (or charge for a negative C-rate).

    Faithful mode reproduces the original program step for step: one linearized solve per
    step, the original write schedule and exit tests, no voltage cutoff. Corrected mode
    solves each step with Newton's method and stops at the voltage cutoffs.
    """
    faithful = p.mode == "faithful"
    if pivot is None:
        pivot = "legacy" if faithful else "partial"
    asm = Assembler(p)
    c = initial_state(p)
    dt = p.dt
    t = 0.0
    mAhg = 0.0
    state = "D" if p.C_rate >= 0 else "C"
    last_write = 0 if faithful else 0.0     # an integer in the original (deviation D-4)
    write_every = p.t_max / p.n_steps / 200
    dc = np.zeros_like(c)
    res = Result()
    n = p.n_steps if max_steps is None else max_steps

    for it in range(1, n + 1):
        row = None
        if it == 1:
            row = _output_row(p, state, t, c, mAhg)
        elif (t - last_write) / 3600 >= write_every:
            row = _output_row(p, state, t, c, mAhg)
            last_write = int(t - dt) if faithful else t
        elif it >= p.n_steps:
            row = _output_row(p, state, t, c, mAhg)
        elif c[-1, P1] >= 99.0 and state == "C":
            res.rows.append(_output_row(p, state, t, c, mAhg))
            res.exit_reason = "end_of_charge"
            break
        elif math.isnan(dc[0, C]):
            res.rows.append(_output_row(p, state, t, c, mAhg))
            res.exit_reason = "nan"
            break
        elif t >= 99.0 * 3600.0:
            res.rows.append(_output_row(p, state, t, c, mAhg))
            res.exit_reason = "max_time"
            break
        if not faithful:
            v = _output_row(p, state, t, c, mAhg)
            if v[2] <= p.V_min or v[2] >= p.V_max:
                res.rows.append(v)
                res.exit_reason = "cutoff_low" if v[2] <= p.V_min else "cutoff_high"
                break
        if row is not None:
            res.rows.append(row)

        # Coulomb counting happens before the solve in the original
        if state == "D":
            mAhg = mAhg + 1000.0 * p.i_specific * dt / 3600.0
        elif state == "C":
            mAhg = mAhg - 1000.0 * p.i_specific * dt / 3600.0

        if faithful:
            A, B, D, G, _ = asm.assemble(c, dt)
            try:
                dc = bandsolver.solve(A, B, D, G, pivot=pivot, backend=backend)
            except (bandsolver.NonFiniteError, bandsolver.SingularBlockError):
                dc = np.full_like(c, np.nan)
            c = c + dc
        else:
            def beyond_cutoff(cn):
                v = _output_row(p, state, t, cn, mAhg)[2]
                return v <= p.V_min or v >= p.V_max
            try:
                c_new, h, stopped = advance(asm, c, dt, backend=backend, pivot=pivot, stop=beyond_cutoff)
            except SolverFailure:
                res.exit_reason = "solver_fail"
                res.steps = it - 1
                break
            dc = c_new - c
            c = c_new
            if stopped:
                sign = 1.0 if state == "D" else -1.0
                mAhg = mAhg - sign * 1000.0 * p.i_specific * (dt - h) / 3600.0
                t = t + h
                v = _output_row(p, state, t, c, mAhg)
                res.rows.append(v)
                res.exit_reason = "cutoff_low" if v[2] <= p.V_min else "cutoff_high"
                res.steps = it
                break
        res.steps = it

        if state == "R":
            dt = dt * 1.0001
        else:
            dt = p.t_max / float(p.n_steps)
        t = t + dt
    else:
        res.exit_reason = res.exit_reason or "max_steps"

    res.final_state = c
    return res


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="Run the LFP model (constant-current discharge).")
    ap.add_argument("--c-rate", type=float, default=1.0)
    ap.add_argument("--mode", choices=("faithful", "corrected"), default="faithful")
    ap.add_argument("--out", default="Time_Voltage.txt")
    args = ap.parse_args(argv)
    p = Params.faithful(C_rate=args.c_rate) if args.mode == "faithful" else Params(C_rate=args.c_rate, mode="corrected")
    r = run(p)
    r.write(args.out)
    print(f"exit: {r.exit_reason} after {r.steps} steps; wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
