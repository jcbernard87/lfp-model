"""Restart files (#23): a run continued from a restart repeats the run it continues."""
import numpy as np
import pytest

from lfp_model.params import Params
from lfp_model.simulate import run

S1, S2 = "cc C=2 t=900", "rest t=300; cc C=-1 Vmax=4.0; cv V=4.0 Imin=0.2"
S2_CV = "cv V=3.35 Imin=0.5"                         # a constant-voltage first step uses the stored current


def _split(tmp_path, s2, **kw):
    """(the rows of S2 in one run of S1; S2, the rows of a run of S2 from a restart written after S1)."""
    whole = run(Params(mode="corrected", steps=f"{S1}; {s2}", **kw)).array
    f = tmp_path / "r.txt"
    run(Params(mode="corrected", steps=S1, restart_write=str(f), **kw))
    part = run(Params(mode="corrected", steps=s2, restart_read=str(f), **kw)).array
    n1 = int(np.sum(whole[:, 7] == 1))                # S1's rows, its end row included
    return whole[n1:], part[1:]                       # a restarted run first repeats the end row


@pytest.mark.parametrize("particle_model", ["uniform", "crystal", "agglomerate"])
@pytest.mark.parametrize("s2", [S2, S2_CV], ids=["cc_rest_cv", "cv_first"])
def test_restart_continues_the_run_exactly(tmp_path, particle_model, s2):
    kw = dict(particle_model=particle_model)
    if particle_model == "agglomerate":
        kw["nj"] = 30
    a, b = _split(tmp_path, s2, **kw)
    assert a.shape == b.shape
    cols = [c for c in range(a.shape[1]) if c != 7]   # every column but the step index (S2's steps renumbered)
    np.testing.assert_array_equal(a[:, cols], b[:, cols])


def test_restart_must_match_the_input(tmp_path):
    f = tmp_path / "r.txt"
    run(Params(mode="corrected", steps="cc C=1 t=60", restart_write=str(f)))
    with pytest.raises(ValueError, match="does not match the input: model"):
        run(Params(mode="corrected", particle_model="crystal", steps="rest t=60", restart_read=str(f)))
    with pytest.raises(ValueError, match="does not match the input: array c"):
        run(Params(mode="corrected", nj=51, steps="rest t=60", restart_read=str(f)))
    with pytest.raises(ValueError, match="mode='corrected'"):
        Params(restart_write=str(f))


@pytest.mark.parametrize("particle_model", ["uniform", "crystal", "agglomerate"])
def test_native_restart_continues_the_run_exactly(fortran_exe, cpp_exe, run_native, tmp_path, particle_model):
    """Fortran and C++: S1 then S2 from a restart equals S1; S2 in one run (byte for byte), and the restart files
    of the two programs and of Python are interchangeable."""
    kw = dict(mode="corrected", particle_model=particle_model, nj=30)
    for exe, tag in ((fortran_exe, "f"), (cpp_exe, "c")):
        whole = run_native(exe, name=f"{tag}w.txt", steps=f"{S1}; {S2}", **kw).read_text().splitlines()
        rf = tmp_path / f"{tag}r.txt"
        run_native(exe, name=f"{tag}1.txt", steps=S1, restart_write=str(rf), **kw)
        part = run_native(exe, name=f"{tag}2.txt", steps=S2, restart_read=str(rf), **kw).read_text().splitlines()
        rows_w = [r.split() for r in whole[2:]]
        n1 = sum(1 for r in rows_w if r[8] == "1")
        a = [r[:8] + r[9:] for r in rows_w[n1:]]
        b = [r.split() for r in part[3:]]
        b = [r[:8] + r[9:] for r in b]
        assert a == b, tag
    # a restart written by one implementation starts the others (Python reads the Fortran program's)
    p = Params(steps=S2, restart_read=str(tmp_path / "fr.txt"), **kw)
    py = run(p).array
    c2 = np.array([[float(x) for x in r.split()[1:]] for r in (tmp_path / "c2.txt").read_text().splitlines()[2:]])
    np.testing.assert_allclose(py[:, :3], c2[:, :3], rtol=0, atol=1e-5)


@pytest.mark.parametrize("which", ["fortran", "cpp"])
def test_native_restart_must_match_the_input(fortran_exe, cpp_exe, run_native, tmp_path, which):
    import subprocess
    exe = fortran_exe if which == "fortran" else cpp_exe
    f = tmp_path / "r.txt"
    run_native(exe, name="a.txt", mode="corrected", steps="cc C=1 t=60", restart_write=str(f))
    for kw, msg in ((dict(particle_model="crystal"), "does not match the input: model"),
                    (dict(nj=51), "does not match the input: array c")):
        with pytest.raises(subprocess.CalledProcessError) as e:
            run_native(exe, name="b.txt", mode="corrected", steps="rest t=60", restart_read=str(f), **kw)
        assert msg in e.value.stderr.decode(), e.value.stderr.decode()
