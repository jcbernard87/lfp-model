from pathlib import Path

import pytest

from lfp_model.namelist import load, parse
from lfp_model.params import Params

ROOT = Path(__file__).resolve().parents[1]


def test_parse_syntax():
    d = parse("""
    ! comment
    &g
      a = 1.5d-3, b = 22  ! trailing
      s = 'x ! not a comment', t = .true.
    /
    """)
    assert d == {"g": {"a": 1.5e-3, "b": 22, "s": "x ! not a comment", "t": True}}


def test_default_input_matches_corrected_defaults():
    p, extra = load(ROOT / "input" / "default.nml")
    assert extra["file"] == "Time_Voltage.txt"
    assert p == Params(mode="corrected")


def test_faithful_example_matches_faithful_defaults():
    p, _ = load(ROOT / "input" / "examples" / "faithful.nml")
    assert p == Params.faithful()


def test_unknown_names_rejected(tmp_path):
    f = tmp_path / "bad.nml"
    f.write_text("&cell\n  nonsense = 1\n/\n")
    with pytest.raises(ValueError, match="nonsense"):
        load(f)
