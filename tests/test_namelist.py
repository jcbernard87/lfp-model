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


# ------------------------------------------------------------------ crystal model inputs

def test_crystal_defaults():
    from lfp_model.params import CRYSTAL_SHAPES
    p = Params(mode="corrected")
    assert (p.particle_model, p.crystal_shape, p.D_c, p.nj_crystal) == ("uniform", "sphere", 8.0e-14, 21)
    assert CRYSTAL_SHAPES == {"slab": 0, "cylinder": 1, "sphere": 2}


@pytest.mark.parametrize("kw, msg", [
    (dict(mode="corrected", particle_model="crystals"), "particle_model"),
    (dict(mode="corrected", particle_model="crystal", crystal_shape="cube"), "crystal_shape"),
    (dict(mode="faithful", particle_model="crystal"), "corrected"),
    (dict(mode="corrected", particle_model="crystal", nj_crystal=3), "nj_crystal"),
])
def test_crystal_inputs_rejected(kw, msg):
    with pytest.raises(ValueError, match=msg):
        Params(**kw)


def test_crystal_keys_read_from_namelist(tmp_path):
    f = tmp_path / "c.nml"
    f.write_text("&cell nj_crystal = 11 /\n&active particle_model = 'crystal', crystal_shape = 'slab', D_c = 1.0d-13 /\n"
                 "&numerics mode = 'corrected' /\n")
    p, _ = load(f)
    assert (p.particle_model, p.crystal_shape, p.D_c, p.nj_crystal) == ("crystal", "slab", 1.0e-13, 11)
