import os
import re
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _exe(env_name, default):
    p = Path(os.environ.get(env_name, ROOT / default))
    if not p.is_absolute():
        p = (Path.cwd() / p).resolve()
    return p if p.is_file() and os.access(p, os.X_OK) else None


@pytest.fixture
def run_native(tmp_path):
    """Run a compiled implementation on the default input with overrides; return the output path."""

    def _run(exe, *, C_rate=1.0, mode="faithful", name="out.txt", steps=None):
        text = (ROOT / "input" / "default.nml").read_text()
        if steps:
            text, n = re.subn(r"steps = '[^']*'", f"steps = '{steps}'", text)
            assert n == 1
        text = re.sub(r"C_rate = [0-9.]+", f"C_rate = {C_rate}", text)
        text = re.sub(r"mode = '\w+'", f"mode = '{mode}'", text)
        text = text.replace("file = 'Time_Voltage.txt'", f"file = '{name}'")
        nml = tmp_path / f"in_{name}.nml"
        nml.write_text(text)
        subprocess.run([str(exe), str(nml)], cwd=tmp_path, check=True, capture_output=True)
        return tmp_path / name

    return _run


@pytest.fixture
def fortran_exe():
    exe = _exe("LFP_FORTRAN_EXE", "build/fortran/lfp_f")
    if exe is None:
        pytest.skip("Fortran program not built (cmake --build build) and LFP_FORTRAN_EXE not set")
    return exe


@pytest.fixture
def cpp_exe():
    exe = _exe("LFP_CPP_EXE", "build/cpp/lfp_cpp")
    if exe is None:
        pytest.skip("C++ program not built (cmake --build build) and LFP_CPP_EXE not set")
    return exe
