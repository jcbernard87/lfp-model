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
    ("&active particle_model = 'crystal', crystal_shape = 'cube' /\n&numerics mode = 'corrected' /\n",
     "crystal_shape must be one of"),
    ("&cell nj_crystal = 3 /\n&active particle_model = 'crystal' /\n&numerics mode = 'corrected' /\n",
     "nj_crystal must be at least 4"),
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


CRYSTAL_CCCV_42 = "cc C=2 Vmin=2.5; rest t=100; cc C=-1 Vmax=4.2; cv V=4.2 Imin=0.5; rest t=10"


@pytest.mark.parametrize("which", ["fortran", "cpp"])
def test_crystal_cccv_to_default_vmax(fortran_exe, cpp_exe, run_native, which):
    """A CC-CV charge held at the default V_max (4.2 V), default D_c: the CV hold starts with the crystal
    surfaces drained to theta ~ 1e-8. Every step completes and every CV row sits on the set voltage."""
    exe = fortran_exe if which == "fortran" else cpp_exe
    out = run_native(exe, mode="corrected", particle_model="crystal", steps=CRYSTAL_CCCV_42, name=f"{which}_cccv.txt")
    _, v = read_tv(out)
    assert int(v[-1, 7]) == 5                                 # the final rest was reached: the protocol completed
    cv = v[v[:, 7] == 4]
    np.testing.assert_allclose(cv[:, 1], 4.2, atol=1e-6)


@pytest.mark.parametrize("C", [2, 5])
def test_crystal_cv_never_accepts_off_setpoint(fortran_exe, cpp_exe, run_native, C):
    """D_c = 1e-17: the crystals cannot sustain the charge, so the CV current search runs into currents
    for which no step converges. A CV state is only accepted on the set voltage; otherwise the run stops
    with a limit reason, the same in both programs."""
    steps = f"cc C={C} Vmin=2.5; rest t=100; cc C=-{C} Vmax=4.2; cv V=4.2 Imin=0.01; rest t=600"
    outs = {}
    for which, exe in (("fortran", fortran_exe), ("cpp", cpp_exe)):
        out = run_native(exe, mode="corrected", particle_model="crystal", D_c=1.0e-17, steps=steps,
                         name=f"{which}_cv{C}.txt")
        _, v = read_tv(out)
        if v[-1, 7] != 5:
            v = v[:-1]                      # a run that stops reports the start-of-step state in its last row
        cv = v[v[:, 7] == 4]
        if len(cv):
            np.testing.assert_allclose(cv[:, 1], 4.2, atol=1e-6)
        outs[which] = out.read_text().splitlines()
    assert [x[:-16] for x in outs["fortran"]] == [x[:-16] for x in outs["cpp"]]


# ------------------------------------------------------------------ driver: cutoffs and failures
@pytest.mark.parametrize("which", ["python", "fortran", "cpp"])
def test_discharge_ignores_its_upper_bound(fortran_exe, cpp_exe, run_native, tmp_path, which):
    """A discharge ends at Vmin only: an upper bound below the starting voltage does not stop it
    (the drivers applied both bounds to every cc step and stopped it at once as cutoff_high)."""
    steps = "cc C=1 Vmax=3.0"
    if which == "python":
        r = run(Params(mode="corrected", steps=steps))
        assert r.exit_reason == "cutoff_low"
        out = tmp_path / "py.txt"
        r.write(out)
    else:
        exe = fortran_exe if which == "fortran" else cpp_exe
        if exe is None:
            pytest.skip(f"{which} not built")
        out = run_native(exe, mode="corrected", steps=steps, name=f"{which}_bound.txt")
    _, v = read_tv(out)
    assert v[-1, 1] == pytest.approx(2.5, abs=1e-4)          # columns after State: time, voltage, ...


@pytest.mark.parametrize("which", ["python", "fortran", "cpp"])
def test_a_physical_limit_keeps_the_progress_of_its_last_step(fortran_exe, cpp_exe, run_native, tmp_path, which):
    """When a step cannot be completed (here the particles fill at 5C), the exit row is the last converged
    sub-step, inside the time step (the drivers reported the state at its start, a whole number of dt)."""
    steps = "cc C=5 Vmin=0.5"
    p = Params(mode="corrected", steps=steps)
    if which == "python":
        r = run(p)
        assert r.exit_reason == "particles_full"
        out = tmp_path / "py.txt"
        r.write(out)
    else:
        exe = fortran_exe if which == "fortran" else cpp_exe
        if exe is None:
            pytest.skip(f"{which} not built")
        out = run_native(exe, mode="corrected", steps=steps, name=f"{which}_limit.txt")
    _, v = read_tv(out)
    t_h = v[-1, 0]                                           # printed in hours to 5 decimals
    grid_h = round(t_h * 3600.0 / p.dt) * p.dt / 3600.0      # the nearest whole number of steps
    assert abs(t_h - grid_h) > 2e-5


@pytest.mark.parametrize("which", ["python", "fortran", "cpp"])
def test_cv_hold_after_a_discharge_proceeds_in_sub_steps(fortran_exe, cpp_exe, run_native, tmp_path, which):
    """A 4.2 V hold right after a 2C discharge cannot be held for whole 10 s steps at first; it proceeds in
    sub-steps and ends on the current limit (without sub-stepping it stopped with particles_empty, #11)."""
    steps = "cc C=2 Vmin=2.5; cv V=4.2 Imin=0.05"
    p = Params(mode="corrected", steps=steps, n_steps=3600)
    if which == "python":
        r = run(p)
        assert r.exit_reason == "end_of_protocol"
        out = tmp_path / "py.txt"
        r.write(out)
    else:
        exe = fortran_exe if which == "fortran" else cpp_exe
        if exe is None:
            pytest.skip(f"{which} not built")
        out = run_native(exe, mode="corrected", steps=steps, n_steps=3600, name=f"{which}_cvsub.txt")
    _, v = read_tv(out)
    assert v[-1, 1] == pytest.approx(4.2, abs=1e-6)                    # voltage held
    assert abs(v[-1, 6]) <= 0.05 * p.i_1C * 1e3 * (1 + 1e-9)            # current [mA/cm2] at the limit
