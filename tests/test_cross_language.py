"""The Fortran (and later C++) programs must agree with the Python package."""
import platform
import subprocess

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
    # byte-identical except the last column, Li_Nernst = (RT/F) ln(c/c_ref): while c relaxes back to
    # c_ref it is round-off (~1e-14 mV), which differs in its leading digit between the compilers
    fl, cl = f_out.read_text().splitlines(), c_out.read_text().splitlines()
    assert len(fl) == len(cl)
    assert [x[:-16] for x in fl] == [x[:-16] for x in cl]
    nf = np.array([float(x.split()[-1]) for x in fl[2:]])
    nc = np.array([float(x.split()[-1]) for x in cl[2:]])
    np.testing.assert_allclose(nf, nc, rtol=1e-4, atol=1e-9)
    r = run(Params(mode="corrected", steps=CYCLE))
    p_out = tmp_path / "py_cycle.txt"
    r.write(p_out)
    sf, vf = read_tv(f_out)
    sp, vp = read_tv(p_out)
    assert sf == sp and vf.shape == vp.shape
    np.testing.assert_allclose(vf, vp, rtol=1e-5, atol=1e-9)


# ------------------------------------------------------------------ crystal model

@pytest.mark.parametrize("shape", ["sphere", "slab"])
def test_fortran_crystal_matches_python(fortran_exe, run_native, tmp_path, shape):
    f_out = run_native(fortran_exe, C_rate=2.0, mode="corrected", particle_model="crystal",
                       crystal_shape=shape, name=f"f_x_{shape}.txt")
    r = run(Params(mode="corrected", C_rate=2.0, particle_model="crystal", crystal_shape=shape))
    p_out = tmp_path / "p.txt"
    r.write(p_out)
    sf, vf = read_tv(f_out)
    sp, vp = read_tv(p_out)
    assert sf == sp and vf.shape == vp.shape
    np.testing.assert_allclose(vf, vp, rtol=1e-5, atol=1e-9)


BAD_CRYSTAL_INPUTS = [
    ("&active particle_model = 'crystals' /\n&numerics mode = 'corrected' /\n", "particle_model must be one of"),
    ("&active particle_model = 'crystal', crystal_shape = 'cube' /\n&numerics mode = 'corrected' /\n", "crystal_shape must be one of"),
    ("&cell nj_crystal = 3 /\n&active particle_model = 'crystal' /\n&numerics mode = 'corrected' /\n", "nj_crystal must be at least 4"),
    ("&active particle_model = 'crystal' /\n&numerics mode = 'faithful' /\n", "needs mode='corrected'"),
]


@pytest.mark.parametrize("text, msg", BAD_CRYSTAL_INPUTS)
def test_fortran_rejects_bad_crystal_input(fortran_exe, tmp_path, text, msg):
    f = tmp_path / "bad.nml"
    f.write_text(text)
    out = subprocess.run([str(fortran_exe), str(f)], cwd=tmp_path, capture_output=True, text=True)
    assert out.returncode != 0 and msg in (out.stdout + out.stderr)


@pytest.mark.parametrize("shape", ["sphere", "cylinder"])
def test_cpp_matches_fortran_crystal(cpp_exe, fortran_exe, run_native, shape):
    kw = dict(C_rate=2.0, mode="corrected", particle_model="crystal", crystal_shape=shape)
    a = run_native(cpp_exe, name=f"c_{shape}.txt", **kw)
    b = run_native(fortran_exe, name=f"f_{shape}.txt", **kw)
    assert a.read_bytes() == b.read_bytes()


CRYSTAL_CYCLE = "cc C=2 Vmin=2.5; rest t=600; cc C=-1 Vmax=4.0; cv V=4.0 Imin=0.05; rest t=600"


def test_crystal_cycle_all_languages(fortran_exe, cpp_exe, run_native, tmp_path):
    kw = dict(mode="corrected", particle_model="crystal", D_c=1.0e-12, steps=CRYSTAL_CYCLE)
    f_out = run_native(fortran_exe, name="fx_cycle.txt", **kw)
    c_out = run_native(cpp_exe, name="cx_cycle.txt", **kw)
    fl, cl = f_out.read_text().splitlines(), c_out.read_text().splitlines()
    assert len(fl) == len(cl) and [x[:-16] for x in fl] == [x[:-16] for x in cl]
    r = run(Params(mode="corrected", particle_model="crystal", D_c=1.0e-12, steps=CRYSTAL_CYCLE))
    p_out = tmp_path / "px.txt"
    r.write(p_out)
    sf, vf = read_tv(f_out)
    sp, vp = read_tv(p_out)
    assert sf == sp and vf.shape == vp.shape
    np.testing.assert_allclose(vf, vp, rtol=1e-5, atol=1e-9)


@pytest.mark.parametrize("text, msg", BAD_CRYSTAL_INPUTS)
def test_cpp_rejects_bad_crystal_input(cpp_exe, tmp_path, text, msg):
    f = tmp_path / "bad.nml"
    f.write_text(text)
    out = subprocess.run([str(cpp_exe), str(f)], cwd=tmp_path, capture_output=True, text=True)
    assert out.returncode != 0 and msg in (out.stdout + out.stderr)
