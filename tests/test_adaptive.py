"""Adaptive time stepping (dt_mode = 'adaptive', #18): backward Euler with the local error estimated from a
linear predictor and controlled against rtol/atol; steps land on the output times."""
import platform

import numpy as np
import pytest

from lfp_model.params import Params
from lfp_model.simulate import run

S = "cc C=2 Vmin=2.5; rest t=600"
_REF = {}


def reference():
    """Fixed steps of 0.25 s; rows every 18 s (a multiple of the step, so on the same times as adaptive rows)."""
    if "r" not in _REF:
        _REF["r"] = run(Params(mode="corrected", steps=S, t_max=3600, n_steps=14400))
    return _REF["r"]


CCCV = "cc C=2 Vmin=2.5; rest t=600; cc C=-1 Vmax=4.0; cv V=4.0 Imin=0.05; rest t=600"


def adaptive(**kw):
    kw.setdefault("steps", S)
    return run(Params(mode="corrected", dt_mode="adaptive", **kw))


def discharge_error(a):
    """Largest voltage difference [V] from the reference at the rows both write during the discharge."""
    ref = reference().array
    common = np.intersect1d(np.round(a[a[:, 7] == 1, 0] * 3600, 6), np.round(ref[ref[:, 7] == 1, 0] * 3600, 6))
    common = common[:-1]                                         # the cutoff rows are located, not on the grid
    va = {round(t * 3600, 6): v for t, v in a[:, :2]}
    vr = {round(t * 3600, 6): v for t, v in ref[:, :2]}
    assert len(common) > 50
    return max(abs(va[t] - vr[t]) for t in common)


def test_fixed_is_the_default():
    assert Params().dt_mode == "fixed"
    with pytest.raises(ValueError, match="dt_mode"):
        Params(dt_mode="variable")


@pytest.mark.parametrize("particle_model", ["uniform", "crystal"])
def test_adaptive_rows_land_on_the_output_times(particle_model):
    """Rows are write_interval apart, as in fixed mode; only the located cutoff and the end of a step are not."""
    a = adaptive(particle_model=particle_model).array
    gap = np.diff(a[:, 0] * 3600)
    assert np.sum(np.abs(gap - 18.0) > 1e-9) <= 2
    assert a[-1, 0] * 3600 == pytest.approx(a[a[:, 7] == 1][-1, 0] * 3600 + 600, abs=1e-9)


def test_adaptive_converges_with_rtol_and_takes_fewer_steps():
    loose, tight = adaptive(rtol=1e-3, atol=1e-3), adaptive(rtol=1e-5, atol=1e-5)
    e_loose, e_tight = discharge_error(loose.array), discharge_error(tight.array)
    assert e_tight < 1e-5 and e_tight < e_loose / 5           # [V]; below ~3e-6 V the reference limits it
    assert loose.steps < tight.steps < reference().steps / 4
    assert loose.exit_reason == tight.exit_reason == reference().exit_reason == "end_of_protocol"


def _table(path):
    lines = path.read_text().splitlines()[2:]
    return [l.split() for l in lines]


@pytest.mark.parametrize("steps", [S, CCCV], ids=["discharge", "cccv"])
@pytest.mark.parametrize("particle_model", ["uniform", "crystal"])
def test_adaptive_native_programs_agree(fortran_exe, cpp_exe, run_native, tmp_path, particle_model, steps):
    """Fortran and C++ take the same adaptive steps (byte-identical output on Apple silicon, as in fixed mode,
    except Li_Nernst; elsewhere to round-off), and agree with Python to the printed precision."""
    kw = dict(mode="corrected", particle_model=particle_model, dt_mode="adaptive", steps=steps)
    f = run_native(fortran_exe, name="f.txt", **kw)
    c = run_native(cpp_exe, name="c.txt", **kw)
    ft, ct = _table(f), _table(c)
    assert len(ft) == len(ct)
    fa = np.array([[float(x) for x in r[1:]] for r in ft])
    ca = np.array([[float(x) for x in r[1:]] for r in ct])
    if platform.system() == "Darwin" and platform.machine() == "arm64":
        assert [r[:9] + r[10:] for r in ft] == [r[:9] + r[10:] for r in ct]     # without Li_Nernst
    else:
        np.testing.assert_allclose(fa[:, :8], ca[:, :8], rtol=1e-9, atol=1e-12)
    py = adaptive(particle_model=particle_model, steps=steps).array
    assert py.shape[0] == fa.shape[0]
    np.testing.assert_allclose(fa[:, :3], py[:, :3], rtol=0, atol=1e-5)        # 5 printed decimals


@pytest.mark.parametrize("which", ["fortran", "cpp"])
def test_native_programs_reject_a_bad_dt_mode(fortran_exe, cpp_exe, run_native, which):
    import subprocess
    exe = fortran_exe if which == "fortran" else cpp_exe
    with pytest.raises(subprocess.CalledProcessError) as e:
        run_native(exe, mode="corrected", dt_mode="variable")
    assert "dt_mode must be one of ('fixed', 'adaptive'), got 'variable'" in e.value.stderr.decode()
