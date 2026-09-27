"""The Fortran (and later C++) programs must agree with the Python package."""
import numpy as np
import pytest

from lfp_model.params import Params
from lfp_model.simulate import run


def read_tv(path):
    lines = open(path).read().splitlines()[2:]
    states = [l.split()[0] for l in lines]
    vals = np.array([[float(x) for x in l.split()[1:]] for l in lines])
    return states, vals


def python_tv(tmp_path, C_rate, mode):
    p = Params.faithful(C_rate=C_rate) if mode == "faithful" else Params(C_rate=C_rate, mode=mode)
    r = run(p)
    out = tmp_path / f"py_{mode}_{C_rate}.txt"
    r.write(out)
    return out


@pytest.mark.parametrize("mode", ["faithful", "corrected"])
@pytest.mark.parametrize("C_rate", [2.0, 0.5])
def test_fortran_matches_python(fortran_exe, run_native, tmp_path, C_rate, mode):
    f_out = run_native(fortran_exe, C_rate=C_rate, mode=mode, name=f"f_{mode}_{C_rate}.txt")
    p_out = python_tv(tmp_path, C_rate, mode)
    if mode == "faithful":
        # both reproduce the original exactly, so their files must be identical
        assert f_out.read_bytes() == p_out.read_bytes()
    else:
        sf, vf = read_tv(f_out)
        sp, vp = read_tv(p_out)
        assert sf == sp and vf.shape == vp.shape
        np.testing.assert_allclose(vf, vp, rtol=1e-5, atol=1e-9, equal_nan=True)


@pytest.mark.parametrize("mode", ["faithful", "corrected"])
@pytest.mark.parametrize("C_rate", [2.0, 0.5])
def test_cpp_matches_python(cpp_exe, run_native, tmp_path, C_rate, mode):
    c_out = run_native(cpp_exe, C_rate=C_rate, mode=mode, name=f"c_{mode}_{C_rate}.txt")
    p_out = python_tv(tmp_path, C_rate, mode)
    if mode == "faithful":
        assert c_out.read_bytes() == p_out.read_bytes()
    else:
        sc, vc = read_tv(c_out)
        sp, vp = read_tv(p_out)
        assert sc == sp and vc.shape == vp.shape
        np.testing.assert_allclose(vc, vp, rtol=1e-5, atol=1e-9, equal_nan=True)


@pytest.mark.parametrize("C_rate", [2.0, 0.5])
def test_cpp_matches_fortran_corrected(cpp_exe, fortran_exe, run_native, C_rate):
    """The two compiled programs implement identical arithmetic, so even corrected runs match exactly."""
    a = run_native(cpp_exe, C_rate=C_rate, mode="corrected", name="c.txt")
    b = run_native(fortran_exe, C_rate=C_rate, mode="corrected", name="f.txt")
    assert a.read_bytes() == b.read_bytes()
