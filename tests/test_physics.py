"""Physics and numerics checks. They need no stored data: everything is computed here.

Corrected mode must pass every check. Faithful mode is checked to still show the documented
defects (docs/deviations.md), so that a change which silently alters faithful behaviour fails.
"""
import numpy as np
import pytest

import bandsolver
from lfp_model import kinetics
from lfp_model.model import Assembler, make_mesh
from lfp_model.params import Params
from lfp_model.simulate import initial_state, newton_step, run

from .reference_residual import faces, residual


def corrected(**kw):
    return Params(mode="corrected", **kw)


def inventories(p, c):
    """Salt in the electrolyte and lithium in the solid, per unit area [mol/cm2]."""
    m = make_mesh(p)
    s, dx = m.s, m.dx
    salt = (p.eps_sep * c[1:s, 0] * dx[1:s]).sum() + (p.eps * c[s + 1:-1, 0] * dx[s + 1:-1]).sum()
    solid = (p.vf_AM * c[s + 1:-1, 3] * dx[s + 1:-1]).sum()
    return salt, solid


@pytest.fixture(scope="module")
def mid_discharge():
    """Corrected-mode state after 600 s at 1C, and the previous state."""
    p = corrected(C_rate=1.0)
    c_prev = run(p, max_steps=599).final_state
    c = run(p, max_steps=600).final_state
    return p, c_prev, c


# ------------------------------------------------------------------ residual and Jacobian

def test_residual_matches_independent_reference(mid_discharge):
    p, c_old, c = mid_discharge
    rng = np.random.default_rng(1)
    cp = c * (1 + 1e-4 * rng.standard_normal(c.shape))       # away from the solution
    asm = Assembler(p)
    _, _, _, G, _ = asm.assemble(cp, 1.0)
    G = G - asm.time_terms(1.0) * (cp - c_old)
    R = residual(p, cp, c_old, 1.0)
    rel = np.abs(G + R).max(axis=0) / np.abs(R).max(axis=0)
    assert rel.max() < 1e-12, rel


def test_faithful_residual_shows_D2_and_D11(mid_discharge):
    """Faithful mode differs from the reference exactly in the cation row (D-2) and current row (D-11)."""
    _, c_old, c = mid_discharge
    p = Params.faithful(C_rate=1.0)
    asm = Assembler(p)
    _, _, _, G, _ = asm.assemble(c, 1.0)
    G = G - asm.time_terms(1.0) * (c - c_old)
    R = residual(p.with_(mode="corrected"), c, c_old, 1.0)
    err = np.abs(G + R).max(axis=0) / np.abs(R).max(axis=0)
    assert err[0] > 0.1 and err[2] > 0.1           # D-2 (cation flux) and D-11 (current)


def _fill(asm, c_old, dt=1.0):
    def fill(c):
        A, B, D, G, _ = asm.assemble(c, dt)
        return A, B, D, G - asm.time_terms(dt) * (c - c_old)
    return fill


def test_jacobian_corrected(mid_discharge):
    p, c_old, c = mid_discharge
    chk = bandsolver.check_jacobian(_fill(Assembler(p), c_old), c)
    assert chk.max_error < 1e-4, chk.worst()


def test_jacobian_faithful_shows_D1(mid_discharge):
    """The Li-foil solid-potential row has the wrong sign in faithful mode (D-1)."""
    _, c_old, c = mid_discharge
    chk = bandsolver.check_jacobian(_fill(Assembler(Params.faithful(C_rate=1.0)), c_old), c)
    name, worst = chk.worst()
    assert (name, worst.node, worst.row, worst.col) == ("B", 0, 1, 1)
    assert worst.user == pytest.approx(-worst.fd, rel=1e-6)


# ------------------------------------------------------------------ conservation

def test_conservation_corrected():
    p = corrected(C_rate=1.0)
    c0 = initial_state(p)
    c = run(p, max_steps=1800).final_state
    s0, so0 = inventories(p, c0)
    s1, so1 = inventories(p, c)
    q = p.i_app * 1800 / p.F
    assert abs(s1 - s0) / s0 < 1e-12                  # salt is conserved
    assert abs((so1 - so0) / q - 1) < 1e-12           # lithium into the solid equals It/F
    N_plus, i2, N_minus = faces(p, c)
    sep = slice(1, make_mesh(p).s - 1)
    np.testing.assert_allclose(i2[sep], p.i_app, rtol=1e-9)   # the separator carries the applied current
    assert np.abs(N_minus[sep]).max() < 1e-6 * p.i_app / p.F  # and no anion flux (quasi-steady)


def test_conservation_faithful_shows_D2():
    p = Params.faithful(C_rate=1.0)
    c0 = initial_state(p)
    c = run(p, max_steps=1800).final_state
    s0, _ = inventories(p, c0)
    s1, _ = inventories(p, c)
    assert -6e-3 < (s1 - s0) / s0 < -4.5e-3            # 0.52 % of the salt is lost at start-up


# ------------------------------------------------------------------ equilibrium

def test_rest_at_equilibrium_stays_put():
    """No current, open-circuit potential everywhere: one step must change nothing."""
    p0 = corrected(C_rate=0.0)
    u = float(kinetics.ocp(p0, np.array([p0.cs_init]))[0])
    p = p0.with_(phi1_init=u)
    asm = Assembler(p)
    c0 = initial_state(p)
    c, _ = newton_step(asm, c0, 10.0)
    assert np.abs(c - c0).max() < 1e-12


# ------------------------------------------------------------------ convergence

def test_time_convergence_first_order():
    """Backward Euler: the Li-face concentration at t = 8 s converges at first order in dt."""
    vals = []
    for dt in (1.0, 0.5, 0.25, 0.125):
        p = corrected(C_rate=2.0, n_steps=int(36000 / dt), newton_tol=1e-13)
        vals.append(run(p, max_steps=int(8 / dt)).final_state[0, 0])
    d = np.diff(vals)
    orders = np.log2(np.abs(d[:-1] / d[1:]))
    assert np.all((orders > 0.9) & (orders < 1.1)), orders


def test_mesh_convergence_second_order():
    """The finite-volume discretization converges at second order in the mesh size."""
    vals = []
    for sep, nj in ((7, 26), (12, 51), (22, 101), (42, 201)):
        p = corrected(C_rate=2.0, sep_node=sep, nj=nj)
        vals.append(run(p, max_steps=900).final_state[-1, 1])
    d = np.diff(vals)
    orders = np.log2(np.abs(d[:-1] / d[1:]))
    assert np.all((orders > 1.8) & (orders < 2.3)), orders


# ------------------------------------------------------------------ end of discharge

@pytest.mark.parametrize("C_rate", [2.0, 1.0])
def test_corrected_ends_at_cutoff(C_rate):
    r = run(corrected(C_rate=C_rate))
    assert r.exit_reason == "cutoff_low"
    v = r.array[:, 1]
    assert v[-1] <= 2.5 and v[-2] > 2.5
    assert 0.79 < r.array[-1, 2] < 0.7977     # electron equivalents: nearly full use of the capacity


# ------------------------------------------------------------------ reference electrode and OCP

def test_ocp_nernst_term_corrected_only():
    """Corrected mode adds (RT/F) ln(c/c_bulk) to the LFP OCP; faithful mode keeps the fit alone."""
    cs = np.array([0.3, 0.5, 0.7]) * kinetics.cs_max(corrected())
    for p in (corrected(), Params.faithful()):
        u0 = kinetics.ocp(p, cs)
        u2 = kinetics.ocp(p, cs, 2.0 * p.c_bulk)
        shift = p.R * p.T / p.F * np.log(2.0) if p.mode == "corrected" else 0.0
        np.testing.assert_allclose(u2 - u0, shift, rtol=1e-12, atol=0.0)


def test_potentials_are_gauge_invariant(mid_discharge):
    """The residual depends only on potential differences, so the foil reference is an exact shift."""
    from lfp_model.simulate import foil_referenced
    p, c_prev, c = mid_discharge
    rng = np.random.default_rng(1)
    x = c * (1.0 + 1.0e-3 * rng.standard_normal(c.shape))    # off the solution, so the residual is not round-off
    r0 = residual(p, x, c_prev, p.dt)
    r1 = residual(p, foil_referenced(p, x, p.i_app), c_prev, p.dt)
    r0[0, 2] = r1[0, 2] = 0.0                                 # the gauge row phi2(foil face) = 0 itself
    np.testing.assert_allclose(r1, r0, rtol=0, atol=1e-9 * np.abs(r0).max())


def test_cell_voltage_against_foil(mid_discharge):
    """V = phi1(collector) - (U_Li + eta_Li): phi2 at the foil face is -(U_Li + eta_Li) on the foil scale."""
    from lfp_model.simulate import cell_voltage, foil_referenced
    p, _, c = mid_discharge
    u_li, eta_li = kinetics.li_foil(p, c[0, 0], p.i_app)
    ref = foil_referenced(p, c, p.i_app)
    assert ref[0, 2] == pytest.approx(-(u_li + eta_li), abs=1e-15)
    assert cell_voltage(p, c, p.i_app) == pytest.approx(ref[-1, 1], abs=1e-15)
    assert eta_li > 0 and c[0, 0] > p.c_bulk and u_li > 0     # discharge: salt builds up at the foil
