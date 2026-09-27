"""Time integration and output (docs/model.md sections 7 and 9, docs/protocol.md).

Faithful mode reproduces the original constant-current discharge: one linearized solve
per step, output rows written before each step, and the original exit tests.
Corrected mode runs a protocol of cc / cv / rest steps, each time step solved with Newton.
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
from .protocol import Step, expand

HEADER = ("State", "Time", "Voltage", "Equivalence", "Anode_Eta", "anode_exchange_c", "Edge_c0")
UNITS = ("CDR", "hours", "Volts", "electron_equivs", "mV", "mA/cm2", "mol/cm3")
HEADER_EXTRA = ("Current", "Step")          # corrected mode only
UNITS_EXTRA = ("mA/cm2", "#")


def _fmt_fixed(v: float) -> str:
    return f"{'NaN':>12}" if math.isnan(v) else f"{v:12.5f}"


def _fmt_sci(v: float) -> str:
    return f"{'NaN':>15}" if math.isnan(v) else f"{v:15.5E}"


def format_header(extended: bool = False) -> str:
    def row(cols):
        # Fortran A<w> right-justifies shorter strings and keeps the leftmost w characters of longer ones
        return (f"{cols[0][:5]:>5} " + " ".join(f"{c[:12]:>12}" for c in cols[1:3]) + " "
                + " ".join(f"{c[:15]:>15}" for c in cols[3:]))
    h, u = (HEADER + HEADER_EXTRA, UNITS + UNITS_EXTRA) if extended else (HEADER, UNITS)
    return row(h) + "\n" + row(u) + "\n"


def format_row(state: str, values) -> str:
    t, v, *rest = values
    cols = [f"{x:15d}" if isinstance(x, int) else _fmt_sci(x) for x in rest]
    return f"{state:>5} {_fmt_fixed(t)} {_fmt_fixed(v)} " + " ".join(cols) + "\n"


@dataclass
class Result:
    rows: list = field(default_factory=list)      # (state, t_h, V, equiv, eta_mV, i0_mA, c_edge[, I_mA, step])
    exit_reason: str = ""
    steps: int = 0
    final_state: Optional[np.ndarray] = None

    @property
    def array(self) -> np.ndarray:
        """Numeric columns: t [h], V, equivalents, eta [mV], i0 [mA/cm2], c_edge (, I [mA/cm2], step)."""
        return np.array([r[1:] for r in self.rows], dtype=float)

    def write(self, path) -> None:
        with open(path, "w") as fh:
            fh.write(format_header(extended=bool(self.rows) and len(self.rows[0]) > 7))
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


def newton_step(asm: Assembler, c_old: np.ndarray, dt: float, I: Optional[float] = None, *,
                backend: str = "fortran", pivot: str = "partial") -> tuple[np.ndarray, int]:
    """One backward-Euler step solved to convergence with Newton's method (corrected mode, fixes D-7)."""
    p = asm.p
    T = asm.time_terms(dt)
    scale = np.array([p.c_bulk, 1.0, 1.0, kinetics.cs_max(p)])
    c = c_old.copy()
    for k in range(1, p.newton_max_iter + 1):
        A, B, D, G, _ = asm.assemble(c, dt, I)
        G = G - T * (c - c_old)
        A, B, D, G = equilibrate(A, B, D, G)
        try:
            dc = bandsolver.solve(A, B, D, G, pivot=pivot, backend=backend)
        except (bandsolver.NonFiniteError, bandsolver.SingularBlockError) as e:
            raise SolverFailure(str(e)) from e
        lam = _bounded_step(c, dc, kinetics.cs_max(p))
        c = c + lam * dc
        if lam == 1.0 and np.max(np.abs(dc) / scale) <= p.newton_tol:
            return c, k
    raise SolverFailure(f"Newton did not converge in {p.newton_max_iter} iterations")


def equilibrate(A, B, D, G):
    """Scale every equation by the largest entry of its row in B.

    The solution is unchanged, but rows of very different size (the solid balance near
    c_s = 0 has entries ~1e14 times the others) no longer trip the solver's relative
    singular-pivot test.
    """
    s = np.abs(B).max(axis=2)
    s[s == 0.0] = 1.0
    return A / s[:, :, None], B / s[:, :, None], D / s[:, :, None], G / s


MAX_DPHI = 0.1   # largest potential change per Newton iteration [V]


def _bounded_step(c, dc, csmax, keep=0.9):
    """Largest step length <= 1 keeping 0 < c, 0 < cs < cs_max (moves at most `keep` of the way to a bound)
    and changing no potential by more than MAX_DPHI (Butler-Volmer exponentials make Newton overshoot)."""
    lam = 1.0
    dphi = float(np.max(np.abs(dc[:, P1:P2 + 1])))
    if dphi > MAX_DPHI:
        lam = MAX_DPHI / dphi
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


EVENT_DV = 1.0e-4    # a cutoff crossing is located to within this voltage [V] ...
EVENT_MIN_DT = 1.0e-9  # ... or this sub-step length [s]


def advance(asm: Assembler, c: np.ndarray, dt: float, I: Optional[float] = None, *, backend: str = "fortran",
            pivot: str = "partial", margin=None, min_dt: float = 1.0e-6):
    """Advance by dt with backward Euler, halving the sub-step on Newton failure.

    `margin(c)` is the distance to the nearest voltage cutoff (negative once crossed). A
    sub-step that crosses by more than EVENT_DV is retried with half the length, so the step
    ends within EVENT_DV of the cutoff. Returns (c, time advanced, stopped).
    """
    t_done, h = 0.0, dt
    while t_done < dt:
        h = min(h, dt - t_done)
        try:
            c_new, _ = newton_step(asm, c, h, I, backend=backend, pivot=pivot)
        except SolverFailure:
            if h / 2 < min_dt:
                raise
            h = h / 2
            continue
        if margin is not None:
            m = margin(c_new)
            if m < 0.0:
                if m < -EVENT_DV and h / 2 >= EVENT_MIN_DT:
                    h = h / 2
                    continue
                return c_new, t_done + h, True
        c, t_done = c_new, t_done + h
    return c, t_done, False


def run(p: Params, *, backend: str = "fortran", pivot: Optional[str] = None, max_steps: Optional[int] = None) -> Result:
    """Run the model: the original discharge in faithful mode, the protocol in corrected mode."""
    if p.mode == "faithful":
        return _run_faithful(p, backend=backend, pivot=pivot or "legacy", max_steps=max_steps)
    return _run_protocol(p, backend=backend, pivot=pivot or "partial", max_steps=max_steps)


def _run_faithful(p: Params, *, backend: str, pivot: str, max_steps: Optional[int]) -> Result:
    """The original program's constant-current discharge, step for step."""
    asm = Assembler(p)
    c = initial_state(p)
    dt = p.dt
    t = 0.0
    mAhg = 0.0
    state = "D" if p.C_rate >= 0 else "C"
    last_write = 0                          # an integer in the original (deviation D-4)
    write_every = p.t_max / p.n_steps / 200
    dc = np.zeros_like(c)
    res = Result()
    n = p.n_steps if max_steps is None else max_steps

    for it in range(1, n + 1):
        if it == 1:
            res.rows.append(_output_row(p, state, t, c, mAhg))
        elif (t - last_write) / 3600 >= write_every:
            res.rows.append(_output_row(p, state, t, c, mAhg))
            last_write = int(t - dt)
        elif it >= p.n_steps:
            res.rows.append(_output_row(p, state, t, c, mAhg))
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

        # Coulomb counting happens before the solve in the original
        if state == "D":
            mAhg = mAhg + 1000.0 * p.i_specific * dt / 3600.0
        elif state == "C":
            mAhg = mAhg - 1000.0 * p.i_specific * dt / 3600.0

        A, B, D, G, _ = asm.assemble(c, dt)
        try:
            dc = bandsolver.solve(A, B, D, G, pivot=pivot, backend=backend)
        except (bandsolver.NonFiniteError, bandsolver.SingularBlockError):
            dc = np.full_like(c, np.nan)
        c = c + dc
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


# ------------------------------------------------------------------ corrected mode: protocols

def li_eta(p: Params, c_edge: float, I: float) -> float:
    """Signed overpotential of the lithium counter electrode (symmetric Butler-Volmer, fixes D-12).

    Negative on discharge (I > 0), so that the cell voltage is phi1(collector) + li_eta.
    """
    alpha = 0.5
    i0 = float(kinetics.li_exchange_current(p, c_edge))
    return -(p.R * p.T / (alpha * p.F)) * math.asinh(I / (2.0 * i0))


def cell_voltage(p: Params, c: np.ndarray, I: float) -> float:
    return float(c[-1, P1]) + li_eta(p, float(c[0, C]), I)


def _row(p: Params, t: float, c: np.ndarray, mAhg: float, I: float, step: int):
    state = "D" if I > 0 else ("C" if I < 0 else "R")
    eta = li_eta(p, float(c[0, C]), I)
    i0 = float(kinetics.li_exchange_current(p, c[0, C]))
    return (state, t / 3600.0, float(c[-1, P1]) + eta, mAhg * p.M * 3.6 / p.F, eta * 1.0e3, i0 * 1.0e3,
            float(c[0, C]), I * 1.0e3, step)


CV_TOL = 1.0e-9      # [V]
CV_MAX_ITER = 50


def cv_step(asm: Assembler, c: np.ndarray, h: float, V_set: float, I_guess: float, *, backend: str, pivot: str):
    """One time step at constant voltage: find the current I with V(I) = V_set.

    f(I) = V(I) - V_set decreases with I. A current the cell cannot carry for the whole step
    (Newton fails) is treated as f = +inf when charging and -inf when discharging, which keeps
    f monotone. The root is bracketed (expanding from the previous current, through I = 0) and
    then found by the Illinois variant of regula falsi, with bisection when a bound is infinite.
    """
    p = asm.p
    states = {}

    def f(I):
        try:
            cn, _ = newton_step(asm, c, h, I, backend=backend, pivot=pivot)
        except SolverFailure:
            return math.inf if I < 0 else -math.inf
        states[I] = cn
        return cell_voltage(p, cn, I) - V_set

    I = I_guess
    fI = f(I)
    if abs(fI) <= CV_TOL:
        return states[I], I
    # bracket: a < b with f(a) > 0 > f(b)
    grow = max(abs(I), 1.0e-2 * p.i_1C)
    a = b = None
    fa = fb = None
    for _ in range(60):
        if fI > 0:
            a, fa = I, fI
            if b is not None:
                break
            I = 0.0 if I < 0 else I + grow
        else:
            b, fb = I, fI
            if a is not None:
                break
            I = 0.0 if I > 0 else I - grow
        grow *= 2.0
        fI = f(I)
        if abs(fI) <= CV_TOL:
            return states[I], I
    if a is None or b is None:
        raise SolverFailure("constant-voltage current could not be bracketed")
    side = 0
    for _ in range(200):
        if math.isfinite(fa) and math.isfinite(fb):
            I = (a * fb - b * fa) / (fb - fa)
            if not (a < I < b):
                I = 0.5 * (a + b)
        else:
            I = 0.5 * (a + b)
        fI = f(I)
        if abs(fI) <= CV_TOL or (b - a) <= 1.0e-14 * p.i_1C:
            if I in states:
                return states[I], I
            raise SolverFailure("constant-voltage step: no feasible current at the set voltage")
        if fI > 0:
            a, fa = I, fI
            if side == 1 and math.isfinite(fb):
                fb *= 0.5
            side = 1
        else:
            b, fb = I, fI
            if side == -1 and math.isfinite(fa):
                fa *= 0.5
            side = -1
    raise SolverFailure("constant-voltage current iteration did not converge")


def _run_protocol(p: Params, *, backend: str, pivot: str, max_steps: Optional[int]) -> Result:
    asm = Assembler(p)
    steps = expand(p)
    c = initial_state(p)
    dt = p.dt
    t, mAhg, n_done = 0.0, 0.0, 0
    res = Result()

    def first_current(st: Step) -> float:
        return st.C * p.i_1C if st.kind == "cc" else 0.0

    I = first_current(steps[0])
    res.rows.append(_row(p, t, c, mAhg, I, 1))
    last_write = t
    reason = ""
    for k, st in enumerate(steps, start=1):
        t_step = 0.0
        if st.kind != "cv":
            I = first_current(st)
        while True:
            if max_steps is not None and n_done >= max_steps:
                res.exit_reason, res.steps, res.final_state = "max_steps", n_done, c
                return res
            h = dt if st.t is None else min(dt, st.t - t_step)
            try:
                if st.kind == "cv":
                    c_new, I = cv_step(asm, c, h, st.V, I, backend=backend, pivot=pivot)
                    h_done = h
                    stopped = st.Imin is not None and abs(I) <= st.Imin * p.i_1C
                    why = "current_limit"
                else:
                    margin = None
                    if st.kind == "cc":
                        def margin(cn, I=I, st=st):
                            v = cell_voltage(p, cn, I)
                            return min(v - st.Vmin, st.Vmax - v)
                    c_new, h_done, stopped = advance(asm, c, h, I, backend=backend, pivot=pivot, margin=margin)
                    why = "cutoff_low" if stopped and cell_voltage(p, c_new, I) <= st.Vmin else "cutoff_high"
            except SolverFailure:
                res.rows.append(_row(p, t, c, mAhg, I, k))
                res.exit_reason, res.steps, res.final_state = "solver_fail", n_done, c
                return res
            mAhg = mAhg + 1000.0 * (I / p.mass_area) * h_done / 3600.0
            c, t, t_step, n_done = c_new, t + h_done, t_step + h_done, n_done + 1
            if not np.all(np.isfinite(c)):
                res.rows.append(_row(p, t, c, mAhg, I, k))
                res.exit_reason, res.steps, res.final_state = "nan", n_done, c
                return res
            if stopped or (st.t is not None and t_step >= st.t * (1.0 - 1.0e-12)):
                res.rows.append(_row(p, t, c, mAhg, I, k))
                last_write = t
                reason = why if stopped else "duration"
                break
            if t - last_write >= p.write_interval:
                res.rows.append(_row(p, t, c, mAhg, I, k))
                last_write = t
            if t >= 99.0 * 3600.0:
                res.rows.append(_row(p, t, c, mAhg, I, k))
                res.exit_reason, res.steps, res.final_state = "max_time", n_done, c
                return res
    res.exit_reason = reason if len(steps) == 1 else "end_of_protocol"
    res.steps, res.final_state = n_done, c
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
