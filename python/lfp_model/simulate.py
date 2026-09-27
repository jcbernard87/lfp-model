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


def run(p: Params, *, backend: str = "fortran", pivot: str = "legacy", max_steps: Optional[int] = None) -> Result:
    """Constant-current discharge, as in the original program."""
    asm = Assembler(p)
    c = initial_state(p)
    dt = p.dt
    t = 0.0
    mAhg = 0.0
    state = "D" if p.C_rate >= 0 else "C"
    last_write = 0            # integer in the original (deviation D-4)
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
        except bandsolver.NonFiniteError:
            dc = np.full_like(c, np.nan)
        except bandsolver.SingularBlockError:
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
