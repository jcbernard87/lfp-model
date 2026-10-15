"""Configurable output (#22): chosen time-series columns, and profiles against x in a long table."""
import subprocess

import numpy as np
import pytest

from lfp_model.params import Params
from lfp_model.simulate import run

STEPS = "cc C=2 t=600; rest t=300"
MODELS = {"uniform": {}, "crystal": {}, "agglomerate": dict(nj=30)}


def _table(path):
    lines = path.read_text().splitlines()
    return lines[:2], np.array([[float(v) for v in ln.split()] for ln in lines[2:]])


@pytest.mark.parametrize("model", sorted(MODELS))
def test_profiles_and_columns_in_every_language(fortran_exe, cpp_exe, run_native, tmp_path, model):
    kw = dict(mode="corrected", particle_model=model, steps=STEPS, profile_every=7, profile_times="100 450",
              columns="Voltage, Current, Th_mean_sep" if model != "uniform" else "Voltage, Current", **MODELS[model])
    files = {}
    for exe, tag in ((fortran_exe, "f"), (cpp_exe, "c")):
        pf = tmp_path / f"{tag}_prof.txt"
        files[tag] = (run_native(exe, name=f"{tag}.txt", profile_file=str(pf), **kw), pf)
    # Fortran and C++: the time series identical; the profiles to round-off (their foil reference holds the
    # Li_Nernst term, whose round-off differs between the compilers, so a last digit can differ)
    assert files["f"][0].read_bytes() == files["c"][0].read_bytes()
    hf, af_ = _table(files["f"][1])
    hc, ac_ = _table(files["c"][1])
    assert hf == hc and af_.shape == ac_.shape
    np.testing.assert_allclose(af_, ac_, rtol=1e-12, atol=1e-15, equal_nan=True)
    pp = tmp_path / "p_prof.txt"
    r = run(Params(profile_file=str(pp), **kw))
    r.write(tmp_path / "p.txt")
    for (fn, pn) in ((files["f"][0], tmp_path / "p.txt"), (files["f"][1], pp)):
        assert fn.read_text().splitlines()[:2] == pn.read_text().splitlines()[:2]
    want = ["State", "Time", "Voltage", "Current"] + (["Th_mean_sep"] if model != "uniform" else [])
    assert (tmp_path / "p.txt").read_text().splitlines()[0].split() == want
    _, af = _table(files["f"][1])
    _, ap = _table(pp)
    assert af.shape == ap.shape
    np.testing.assert_allclose(af, ap, rtol=1e-6, atol=1e-12, equal_nan=True)


def test_profile_schedule_and_values(tmp_path):
    """Profiles with every 7th row, at the first row at or after 100 s and 450 s, and with the last row; theta in
    the profile is the time series' Th_mean_sep at the first crystal; c starts at c_bulk everywhere."""
    pp = tmp_path / "prof.txt"
    p = Params(mode="corrected", particle_model="crystal", steps=STEPS, profile_every=7, profile_times="450, 100",
               profile_file=str(pp))
    r = run(p)
    t_rows = r.array[:, 0] * 3600.0
    due = {i for i in range(len(t_rows)) if i % 7 == 0}
    for tp in (100.0, 450.0):
        due.add(int(np.argmax(t_rows >= tp)))
    due.add(len(t_rows) - 1)
    assert [round(pr[0], 6) for pr in r.profiles] == [round(t_rows[i], 6) for i in sorted(due)]
    _, a = _table(pp)
    assert a.shape[0] == len(due) * p.nj
    first = a[:p.nj]
    assert np.all(first[:, 2] == p.c_bulk)
    s = p.sep_node - 1
    th_mean_sep = r.array[:, 10]
    for k, i in enumerate(sorted(due)):
        prof = a[k * p.nj:(k + 1) * p.nj]
        assert prof[s + 1, 5] == pytest.approx(th_mean_sep[i], rel=1e-6)
        assert np.isnan(prof[0, 3]) and np.isnan(prof[0, 5])          # no solid in the separator


def test_default_output_is_unchanged(tmp_path):
    """With no columns and no profile file the output is exactly as before (the original column set)."""
    r = run(Params(mode="corrected", steps="cc C=1 t=60"))
    r.write(tmp_path / "o.txt")
    head = (tmp_path / "o.txt").read_text().splitlines()[0].split()
    assert head == ["State", "Time", "Voltage", "Equivalence", "Anode_Eta", "anode_exchange_", "Edge_c0", "Current",
                    "Step", "Li_Nernst"]
    assert r.profiles == []


@pytest.mark.parametrize("which", ["python", "fortran", "cpp"])
@pytest.mark.parametrize("kw,msg", [(dict(columns="Voltage, Vltage"), "unknown output column 'Vltage'"),
                                    (dict(profile_file="x.txt", profile_vars="c, theta_surf"),
                                     "unknown profile variable 'theta_surf'")])
def test_unknown_names_are_rejected(fortran_exe, cpp_exe, run_native, tmp_path, which, kw, msg):
    kw = dict(kw)
    if "profile_file" in kw:
        kw["profile_file"] = str(tmp_path / kw["profile_file"])
    if which == "python":
        with pytest.raises(ValueError, match=msg):
            r = run(Params(mode="corrected", steps="cc C=1 t=60", **kw))
            r.write(tmp_path / "o.txt")
        return
    exe = fortran_exe if which == "fortran" else cpp_exe
    with pytest.raises(subprocess.CalledProcessError) as e:
        run_native(exe, mode="corrected", steps="cc C=1 t=60", **kw)
    assert msg in e.value.stderr.decode(), e.value.stderr.decode()
