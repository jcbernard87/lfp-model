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
