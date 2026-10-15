"""The programs' built-in defaults are the documented ones: a minimal input gives the same run as the full
input file, in Fortran, C++ and Python."""
import re
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SHORT = "cc C=1 t=600; rest t=300"


def _full(model, name):
    text = (ROOT / "input" / "default.nml").read_text(encoding="utf-8")
    text = re.sub(r"particle_model = '\w+'", f"particle_model = '{model}'", text)
    text = re.sub(r"steps = '[^']*'", f"steps = '{SHORT}'", text)
    return text.replace("file = 'Time_Voltage.txt'", f"file = '{name}'")


def _minimal(model, name):
    return (f"&active\n  particle_model = '{model}'\n/\n&numerics\n  mode = 'corrected'\n/\n"
            f"&protocol\n  steps = '{SHORT}'\n/\n&output\n  file = '{name}'\n/\n")


@pytest.mark.parametrize("model", ["uniform", "crystal"])
def test_programs_default_to_the_input_file(model, fortran_exe, cpp_exe, tmp_path):
    """A minimal input (particle model, mode, protocol, output file) gives the same run as default.nml."""
    (tmp_path / "full.nml").write_text(_full(model, "full.txt"))
    (tmp_path / "min.nml").write_text(_minimal(model, "min.txt"))
    for exe in (fortran_exe, cpp_exe):
        for nml in ("full.nml", "min.nml"):
            subprocess.run([str(exe), nml], cwd=tmp_path, check=True, capture_output=True)
        assert (tmp_path / "full.txt").read_bytes() == (tmp_path / "min.txt").read_bytes(), exe


@pytest.mark.parametrize("model", ["uniform", "crystal"])
def test_python_defaults_to_the_input_file(model, tmp_path):
    from lfp_model.namelist import load
    (tmp_path / "full.nml").write_text(_full(model, "full.txt"))
    (tmp_path / "min.nml").write_text(_minimal(model, "min.txt"))
    assert load(tmp_path / "full.nml")[0] == load(tmp_path / "min.nml")[0]


BAD = [
    ("&cellz\n  nj = 30\n/\n", "unknown namelist group &cellz"),
    ("&numerics\n  nj = 30\n/\n", ("nj", "&numerics")),                       # a name in the wrong group
    ("&cell\n  nj = 30\n/\n&cell\n  nj = 40\n/\n", "namelist group &cell appears twice"),
]


@pytest.mark.parametrize("text,msg", BAD)
def test_programs_reject_bad_groups(text, msg, fortran_exe, cpp_exe, tmp_path):
    """An unknown group, a name in the wrong group and a repeated group stop every implementation with a
    message naming the problem. The Fortran and C++ programs skipped unknown groups silently, the C++ program
    accepted a name in any group, and a repeated group was read three ways (Fortran first, C++ merged, Python
    last)."""
    from lfp_model.namelist import load
    nml = tmp_path / "bad.nml"
    nml.write_text(text + "&output\n  file = 'out.txt'\n/\n")
    msgs = msg if isinstance(msg, tuple) else (msg,)
    for exe in (fortran_exe, cpp_exe):
        r = subprocess.run([str(exe), str(nml)], cwd=tmp_path, capture_output=True, text=True)
        assert r.returncode != 0 and all(m in r.stdout + r.stderr for m in msgs), (exe, r.stdout + r.stderr)
    with pytest.raises(ValueError) as err:
        load(nml)
    assert all(m in str(err.value) for m in msgs), str(err.value)


def test_a_quoted_path_with_slashes_is_one_value(fortran_exe, cpp_exe, tmp_path):
    """A '/' inside a quoted value (an output path with directories) does not end the namelist group: the C++
    and Python readers ended the group there, so the output file name was lost."""
    from lfp_model.namelist import load
    out = tmp_path / "sub" / "out.txt"
    out.parent.mkdir()
    nml = tmp_path / "in.nml"
    nml.write_text(f"&numerics\n  mode = 'corrected'\n/\n&protocol\n  steps = 'cc C=1 t=60'\n/\n"
                   f"&output\n  file = '{out}', write_interval = 30.0\n/\n")
    for exe in (fortran_exe, cpp_exe):
        out.unlink(missing_ok=True)
        r = subprocess.run([str(exe), str(nml)], cwd=tmp_path, capture_output=True, text=True)
        assert r.returncode == 0 and out.exists(), (exe, r.stdout + r.stderr)
    p, extra = load(nml)
    assert extra["file"] == str(out) and p.write_interval == 30.0


def test_fixture_values_reach_the_input(tmp_path):
    """The native-input fixture sets every value it is given, also names that share a line with others in
    default.nml (an earlier version left those at the default: the line set them again after the new entry)."""
    from lfp_model.namelist import load
    from .conftest import native_input
    over = dict(nj=31, sep_node=12, eps=0.45, n_steps=7200, D=1.5e-6, write_interval=30.0, R_p=1.5e-5,
                particle_model="agglomerate", nja=9, R_agg=4.0e-4, dt_mode="adaptive", rtol=1e-5)
    nml = tmp_path / "in.nml"
    nml.write_text(native_input(mode="corrected", C_rate=2.0, steps="cc C=2 t=60", **over))
    p, _ = load(nml)
    for k, v in over.items():
        assert getattr(p, k) == v, k
    assert (p.mode, p.C_rate, p.steps) == ("corrected", 2.0, "cc C=2 t=60")
    with pytest.raises(KeyError):
        native_input(not_a_name=1)
