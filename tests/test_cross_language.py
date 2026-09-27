"""The Fortran (and later C++) programs must agree with the Python package."""
import platform

import numpy as np
import pytest

from lfp_model.params import Params
from lfp_model.simulate import run


def read_tv(path):
    lines = open(path).read().splitlines()[2:]
    states = [l.split()[0] for l in lines]
    vals = np.array([[float(x) for x in l.split()[1:]] for l in lines])
    return states, vals


EXACT_PLATFORM = platform.system() == "Darwin" and platform.machine() == "arm64"


def assert_faithful_match(native, py):
    """Faithful outputs are byte-identical on the reference platform (macOS arm64, where numpy and
    the compiled programs share libm). Elsewhere numpy's vectorized exp/atan may differ from libm
    in the last bit, so the printed values are compared instead."""
    if EXACT_PLATFORM:
        assert native.read_bytes() == py.read_bytes()
        return
    sn, vn = read_tv(native)
    sp, vp = read_tv(py)
    assert sn == sp and vn.shape == vp.shape
    np.testing.assert_allclose(vn, vp, rtol=2e-5, atol=1e-9, equal_nan=True)


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
        assert_faithful_match(f_out, p_out)
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
        assert_faithful_match(c_out, p_out)
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


CYCLE = "cc C=2 Vmin=2.5; rest t=600; cc C=-1 Vmax=4.0; cv V=4.0 Imin=0.05; rest t=600"


def test_protocol_cycle_all_languages(fortran_exe, cpp_exe, run_native, tmp_path):
    """A full discharge / rest / charge / CV / rest cycle agrees across the three implementations."""
    f_out = run_native(fortran_exe, mode="corrected", name="f_cycle.txt", steps=CYCLE)
    c_out = run_native(cpp_exe, mode="corrected", name="c_cycle.txt", steps=CYCLE)
    assert f_out.read_bytes() == c_out.read_bytes()
    r = run(Params(mode="corrected", steps=CYCLE))
    p_out = tmp_path / "py_cycle.txt"
    r.write(p_out)
    sf, vf = read_tv(f_out)
    sp, vp = read_tv(p_out)
    assert sf == sp and vf.shape == vp.shape
    np.testing.assert_allclose(vf, vp, rtol=1e-5, atol=1e-9)
