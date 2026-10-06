"""Agglomerate model (particle_model = 'agglomerate', docs/model.md section 13), ported from nmc111-model's
corrected agglomerate model: conservation, the condensed Newton step, limits, convergence."""
import math

import numpy as np
import pytest

from lfp_model.agglomerate import AggState, AgglomerateModel
from lfp_model.logcore import P1, P2, S, U, sigmoid
from lfp_model.params import Params
from lfp_model.simulate import run


def params(**kw):
    dt = kw.pop("dt", 5.0)
    base = dict(mode="corrected", particle_model="agglomerate", nj=30, nja=10, t_max=36000.0,
                n_steps=int(round(36000.0 / dt)))
    base.update(kw)
    return Params(**base)


def salt(m, p, st):
    """Salt in the macro-pores and in the agglomerate pores [mol/cm2]."""
    me = m.mesh
    s = me.s
    ce, ca = m.conc(st), p.c_bulk * np.exp(st.ca[..., U])
    macro = (p.eps_sep * ce[1:s] * me.dx[1:s]).sum() + (p.eps * ce[s + 1:-1] * me.dx[s + 1:-1]).sum()
    v_one = 4.0 / 3.0 * math.pi * p.R_agg ** 3
    micro = ((p.eps_agg * ca * m.dV).sum(axis=1) / v_one * m.v_agg * me.dx[m.nodes]).sum()
    return macro + micro


def lithium_in_crystals(m, p, st):
    v_one = 4.0 / 3.0 * math.pi * p.R_agg ** 3
    per_agg = ((m.cs(st) - p.cs_init) * m.dV).sum(axis=1) * (1.0 - p.eps_agg) / v_one * m.v_agg
    return (per_agg * m.mesh.dx[m.nodes]).sum()


@pytest.fixture(scope="module")
def discharge():
    p = params(steps="cc C=1 Vmin=2.5")
    return p, run(p)


def test_discharge_ends_at_cutoff(discharge):
    p, r = discharge
    assert isinstance(r.final_state, AggState)            # the agglomerate model actually ran
    assert r.exit_reason == "cutoff_low"
    v = r.array[:, 1]
    assert v[-1] <= 2.5 < v[-2] and abs(v[-1] - 2.5) < 1e-4


def test_lithium_conserved(discharge):
    """Lithium in the crystals equals the charge passed (1C throughout), to round-off."""
    p, r = discharge
    m = AgglomerateModel(p)
    passed = r.array[-1, 0] * 3600.0 * p.i_1C / p.F
    assert lithium_in_crystals(m, p, r.final_state) / passed == pytest.approx(1.0, abs=1e-12)


def test_salt_conserved(discharge):
    """Salt in the macro-pores plus the agglomerate pores is constant (Li+ in at the foil = Li+ into the crystals)."""
    p, r = discharge
    m = AgglomerateModel(p)
    s0 = salt(m, p, m.initial_state())
    assert abs(salt(m, p, r.final_state) - s0) / s0 < 1e-13


def test_condensed_step_solves_the_full_system():
    """After a condensed Newton step, both scales' uncondensed residuals vanish (each residual column compared
    with its size at a randomly perturbed state)."""
    p = params()
    m = AgglomerateModel(p)
    old = m.initial_state()
    st = m.newton_step(old, p.dt, p.i_1C)
    rng = np.random.default_rng(7)
    x = AggState(st.c + 1e-3 * rng.standard_normal(st.c.shape), st.ca + 1e-3 * rng.standard_normal(st.ca.shape))
    Re, Ra = m.residuals(st, old, p.dt, p.i_1C)
    Rex, Rax = m.residuals(x, old, p.dt, p.i_1C)
    for R, Rx, cols in ((Re, Rex, (U, P1, P2)), (Ra, Rax, (U, P1, P2, S))):
        for col in cols:
            assert np.abs(R[..., col]).max() <= 1e-9 * np.abs(Rx[..., col]).max()


def test_newton_converges_quadratically():
    """From a perturbed iterate, the scaled update falls quadratically: the condensed Jacobian is exact."""
    p = params(newton_tol=1e-14)
    m = AgglomerateModel(p)
    old = m.initial_state()
    sol = m.newton_step(old, p.dt, p.i_1C)
    rng = np.random.default_rng(3)
    x = AggState(sol.c.copy(), sol.ca.copy())
    x.c[:, P1:P2 + 1] += 1e-3 * rng.standard_normal((p.nj, 2))
    x.ca[..., P1:P2 + 1] += 1e-3 * rng.standard_normal(x.ca[..., P1:P2 + 1].shape)
    x.ca[..., U] += 1e-2 * rng.standard_normal(x.ca[..., U].shape)
    hist = []
    m.newton_step(old, p.dt, p.i_1C, start=x, history=hist)
    assert len(hist) <= 6
    for a, b in zip(hist, hist[1:]):
        if a > 1e-7:
            assert b < 50.0 * a ** 2


def test_charge_neutral_agglomerates(discharge):
    """Each agglomerate takes in as much electronic as it gives ionic current: i1,in + i2,in = 0."""
    p, r = discharge
    q, _ = AgglomerateModel(p).surface_flux(r.final_state.ca)
    assert np.abs(q[:, 1] + q[:, 2]).max() <= 1e-7 * np.abs(q[:, 1]).max()   # Newton tolerance level


def test_rest_at_equilibrium_stays_put():
    p = params()
    m = AgglomerateModel(p)
    st = m.initial_state()
    u = float(m.kin.ocp(0.0, st.ca[0, 0, S]))
    st.c[:, P1] = u
    st.ca[..., P1] = u
    new = m.newton_step(st, 10.0, 0.0)
    assert np.abs(new.c - st.c).max() < 1e-12 and np.abs(new.ca - st.ca).max() < 1e-12


def test_uniform_limit():
    """With fast transport and conduction inside the agglomerates and a vanishing pore volume (whose electrolyte
    the uniform model does not have), the agglomerate model becomes the uniform-particle model: max |dV| is
    2.1e-9 V at eps_agg = 1e-6 and 4.9e-12 V at 1e-8. Fast transport alone makes the crystals lithiate evenly."""
    S_ = "cc C=5 t=300"
    u = run(Params(mode="corrected", steps=S_)).array
    p = Params(mode="corrected", particle_model="agglomerate", steps=S_, eps_agg=1e-6, tau_agg=1e-12, sigma_agg=1e3)
    a = run(p)
    assert np.abs(a.array[:, 1] - u[:, 1]).max() < 1e-8
    fast = Params(mode="corrected", particle_model="agglomerate", steps=S_, tau_agg=1e-4, sigma_agg=1e3)
    th = AgglomerateModel(fast).cs(run(fast).final_state)
    assert float(((th.max(axis=1) - th.min(axis=1)) / th.mean(axis=1)).max()) < 1e-3


def test_slow_agglomerate_transport_costs_voltage():
    """A large agglomerate, tortuous pores or a poorly conducting agglomerate solid lowers the discharge voltage."""
    S_ = "cc C=5 t=300"
    v = lambda **kw: run(Params(mode="corrected", particle_model="agglomerate", steps=S_, **kw)).array[-1, 1]
    base = v()
    assert v(R_agg=2e-3) < base and v(tau_agg=50.0) < base and v(sigma_agg=1e-6) < base - 0.02


CYCLE = "cc C=1 Vmin=2.5; rest t=600; cc C=-1 Vmax=4.0; cv V=4.0 Imin=0.05; rest t=600"


def test_protocol_cycle_returns_the_lithium():
    p = params(steps=CYCLE)
    r = run(p)
    assert r.exit_reason == "end_of_protocol"
    li = lithium_in_crystals(AgglomerateModel(p), p, r.final_state)
    net = r.array[-1, 2] * p.mass_area / p.M        # electron equivalents x moles of active material
    assert li == pytest.approx(net, rel=1e-9, abs=1e-15)


def _v600(**kw):
    """Cell voltage after 600 s at 2C."""
    base = dict(nj=30, nja=10, dt=10.0, steps="cc C=2 t=600")
    base.update(kw)
    return run(params(**base)).array[-1, 1]


def _orders(values):
    e = np.abs(np.diff(values))
    return [math.log2(e[i] / e[i + 1]) for i in range(len(e) - 1)]


def test_time_convergence_first_order():
    """Backward Euler: halving the step halves the error, in u = ln(c/c_bulk) at the foil after 8 s at 2C (orders
    0.96, 0.98, 0.99 measured). Later in the discharge the cell is quasi-steady on the flat LiFePO4 plateau and
    the time error, ~1e-8, is too small to measure an order."""
    def u_foil(dt):
        return run(params(dt=dt, steps="cc C=2 t=8", newton_tol=1e-13)).final_state.c[0, U]
    orders = _orders([u_foil(d) for d in (1.0, 0.5, 0.25, 0.125)])
    assert all(0.9 < o < 1.1 for o in orders), orders


def test_radial_mesh_convergence_second_order():
    """Halving the agglomerate's radial spacing (4, 8, 16, 32 volumes) quarters the error."""
    orders = _orders([_v600(nja=n, sigma_agg=1e-5) for n in (6, 10, 18, 34)])
    assert all(1.8 < o < 2.2 for o in orders), orders


def test_electrode_mesh_convergence_second_order():
    """Halving the cathode spacing (8, 16, 32, 64 volumes) quarters the error."""
    orders = _orders([_v600(nj=22 + 1 + n) for n in (8, 16, 32, 64)])
    assert all(1.8 < o < 2.2 for o in orders), orders


@pytest.mark.parametrize("C_rate,tau_agg", [(1.0, 22400.0), (3.0, 2240.0)])
def test_exhausted_agglomerate_cores_run_to_the_cutoff(C_rate, tau_agg):
    """With slow transport inside the agglomerates their cores run out of salt; the discharge still ends smoothly
    at the voltage cutoff, with lithium and salt conserved."""
    p = params(dt=10.0, tau_agg=tau_agg, steps=f"cc C={C_rate} Vmin=2.5")
    r = run(p)
    assert r.exit_reason == "cutoff_low"
    m = AgglomerateModel(p)
    st = r.final_state
    assert (p.c_bulk * np.exp(st.ca[..., U])).min() < 1e-3 * p.c_bulk          # the cores are drained
    passed = r.array[-1, 0] * 3600.0 * C_rate * p.i_1C / p.F
    assert lithium_in_crystals(m, p, st) / passed == pytest.approx(1.0, abs=1e-12)
    s0 = salt(m, p, m.initial_state())
    assert abs(salt(m, p, st) - s0) / s0 < 1e-12


def test_charge_empties_the_crystals_smoothly():
    p = params(dt=10.0, steps="cc C=1 Vmin=2.5; rest t=600; cc C=-1 Vmax=4.2; cv V=4.2 Imin=0.02; rest t=600")
    r = run(p)
    assert r.exit_reason == "end_of_protocol"
    m = AgglomerateModel(p)
    assert (m.cs(r.final_state) / m.cs_max).max() < 1e-3        # nearly empty, reached without clipping


def test_rest_relaxes_monotonically():
    """After a 2C partial discharge (with a poorly conducting agglomerate solid, so the agglomerates develop
    gradients) the voltage moves monotonically during a rest and the pore electrolyte relaxes. The spread of theta
    does not: on the flat LiFePO4 plateau crystals at different theta are nearly at the same potential (unlike
    NMC111, whose agglomerates relax to uniform theta)."""
    def spreads(st):
        th, u = sigmoid(st.ca[..., S]), st.ca[..., U]
        return float(np.max(th.max(axis=1) - th.min(axis=1))), float(np.max(u.max(axis=1) - u.min(axis=1)))

    before = spreads(run(params(steps="cc C=2 t=900", sigma_agg=1e-5)).final_state)
    r = run(params(steps="cc C=2 t=900; rest t=1800", sigma_agg=1e-5))
    v = np.array([row[2] for row in r.rows if row[8] == 2])
    dv = np.diff(v)
    assert len(v) > 10 and (np.all(dv >= -1e-9) or np.all(dv <= 1e-9)), dv
    after = spreads(r.final_state)
    assert before[1] > 1e-4 and after[1] < 1e-3 * before[1], (before, after)


def test_agglomerate_columns_in_the_output(tmp_path):
    """The agglomerate model writes the surface-crystal and mean theta of the agglomerates next to the separator
    and the collector, under the crystal model's column names."""
    r = run(params(steps="cc C=2 t=600"))
    out = tmp_path / "o.txt"
    r.write(out)
    assert out.read_text().splitlines()[0].split()[-4:] == ["Th_surf_sep", "Th_mean_sep", "Th_surf_col", "Th_mean_col"]
    ts_sep, tm_sep, ts_col, tm_col = r.array[-1, -4:]
    assert 0 < tm_sep <= ts_sep * (1 + 1e-12) and 0 < tm_col <= ts_col * (1 + 1e-12)


# ------------------------------------------------------------------ the Fortran and C++ programs
def _tokens(path):
    return [ln.split() for ln in path.read_text().splitlines()]


@pytest.mark.parametrize("steps,extra", [("cc C=2 Vmin=2.5; rest t=600", {}), (CYCLE, {}),
                                         ("cc C=2 Vmin=2.5; rest t=600", {"dt_mode": "adaptive"})],
                         ids=["discharge", "cycle", "adaptive"])
def test_agglomerate_native_programs_agree(fortran_exe, cpp_exe, run_native, tmp_path, steps, extra):
    """The Fortran and C++ programs write the same files (Li_Nernst, token 9, to round-off), and Python agrees to
    the printed precision, the theta columns included."""
    kw = dict(mode="corrected", particle_model="agglomerate", nj=30, steps=steps, **extra)
    f = run_native(fortran_exe, name="f.txt", **kw)
    c = run_native(cpp_exe, name="c.txt", **kw)
    ft, ct = _tokens(f), _tokens(c)
    assert len(ft) == len(ct) and [r[:9] + r[10:] for r in ft] == [r[:9] + r[10:] for r in ct]
    r = run(Params(**kw))
    py = tmp_path / "p.txt"
    r.write(py)
    pt = _tokens(py)
    assert pt[:2] == ft[:2] and len(pt) == len(ft)
    fa = np.array([[float(x) for x in row[1:]] for row in ft[2:]])
    pa = np.array([[float(x) for x in row[1:]] for row in pt[2:]])
    np.testing.assert_allclose(fa[:, :3], pa[:, :3], rtol=0, atol=1e-5)          # 5 printed decimals
    np.testing.assert_allclose(fa[:, 9:], pa[:, 9:], rtol=1e-5, atol=1e-9)


@pytest.mark.parametrize("extra,msg", [
    (dict(nja=3), "nja must be at least 4, got 3"),
    (dict(eps_agg=0.3), "porosity + agglomerate volume fraction exceed 1"),
    (dict(mode="faithful"), "particle_model='agglomerate' needs mode='corrected'"),
])
def test_bad_agglomerate_inputs_are_rejected(fortran_exe, cpp_exe, run_native, extra, msg):
    import subprocess
    kw = dict(mode="corrected", particle_model="agglomerate")
    kw.update(extra)
    for exe in (fortran_exe, cpp_exe):
        with pytest.raises(subprocess.CalledProcessError) as e:
            run_native(exe, **kw)
        assert msg in e.value.stderr.decode(), (exe, e.value.stderr.decode())
    with pytest.raises(ValueError, match=msg.replace("+", r"\+")):
        Params(**kw)
