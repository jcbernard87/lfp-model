"""Command line: python -m lfp_model [input.nml]"""
import sys
from pathlib import Path

from .namelist import load
from .simulate import run


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    path = Path(argv[0]) if argv else Path(__file__).resolve().parents[2] / "input" / "default.nml"
    p, extra = load(path)
    out = extra.get("file", "Time_Voltage.txt")
    r = run(p)
    r.write(out)
    print(f"{p.mode} run, C-rate {p.C_rate:g}: exit '{r.exit_reason}' after {r.steps} steps; wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
