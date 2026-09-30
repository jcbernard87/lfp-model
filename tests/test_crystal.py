import math
import numpy as np
import pytest

import bandsolver
from lfp_model.crystal import CrystalMesh, diffusion_rows
from lfp_model.logcore import sigmoid, logit


def _run_flux(k, nc, J=0.05, th0=0.2, dt=0.05, t_end=2.0):
    """Unit crystal (R = D = cs_max = 1), constant inward flux J; Newton in s each step."""
    cm = CrystalMesh(1.0, nc, k)
    xc = np.full((1, nc), logit(th0))
    for _ in range(round(t_end / dt)):
        old = xc.copy()
        for _ in range(30):
            R, B, Dn, An = diffusion_rows(cm, 1.0, 1.0, xc, old, dt)
            R[:, -1] += J                       # surface balance R_surf - i/F; an inward flux J is i/F = -J
            A3 = An[0][:, None, None]; B3 = B[0][:, None, None]; D3 = Dn[0][:, None, None]
            dx = bandsolver.solve(A3, B3, D3, -R[0][:, None])[:, 0]
            xc = xc + dx[None, :]
            if np.abs(dx).max() < 1e-13:
                break
    th = sigmoid(xc[0])
    exact = th0 + (k + 1) * J * t_end + J / 2 * (cm.r ** 2 - (k + 1) / (k + 3))
    return cm, th, exact


@pytest.mark.parametrize("k", [0, 1, 2])
def test_constant_flux_matches_analytic(k):
    cm, th, exact = _run_flux(k, 41)
    inner = slice(1, None)                      # node 0 copies node 1 (not at r = 0)
    assert np.abs(th[inner] - exact[inner]).max() < 2e-4
    mean = (cm.V * th).sum() / cm.Vtot
    assert mean == pytest.approx(0.2 + (k + 1) * 0.05 * 2.0, abs=1e-12)    # lithium balance, exact


@pytest.mark.parametrize("k", [0, 2])
def test_crystal_mesh_second_order(k):
    errs = []
    for nc in (12, 22, 42):
        _, th, exact = _run_flux(k, nc)
        errs.append(abs(th[-1] - exact[-1]))    # surface value
    orders = np.log2(np.array(errs[:-1]) / np.array(errs[1:]))   # interior cells 10, 20, 40: h halves
    assert np.all(orders > 1.8), orders


def test_diffusion_rows_jacobian():
    rng = np.random.default_rng(3)
    cm = CrystalMesh(2e-5, 9, 2)
    xc = rng.normal(0.0, 2.0, (2, 9)); old = xc + rng.normal(0.0, 0.1, xc.shape)
    R, B, Dn, An = diffusion_rows(cm, 8e-14, 0.0228, xc, old, 1.0)
    h = 1e-6
    for j in range(9):
        e = np.zeros_like(xc); e[:, j] = h
        dR = (diffusion_rows(cm, 8e-14, 0.0228, xc + e, old, 1.0)[0] - diffusion_rows(cm, 8e-14, 0.0228, xc - e, old, 1.0)[0]) / (2 * h)
        want = np.zeros_like(xc)
        want[:, j] = B[:, j]
        if j > 0: want[:, j - 1] = Dn[:, j - 1]
        if j < 8: want[:, j + 1] = An[:, j + 1]
        np.testing.assert_allclose(dR, want, rtol=1e-6, atol=1e-6 * np.abs(want).max())


from lfp_model.params import Params
from lfp_model.crystal import CrystalModel, CrystalState


def crystal(**kw):
    return Params(mode="corrected", particle_model="crystal", **kw)


def test_crystal_area_by_shape():
    for shape, k in (("slab", 0), ("cylinder", 1), ("sphere", 2)):
        p = crystal(crystal_shape=shape)
        assert CrystalModel(p).a == pytest.approx((k + 1) * p.vf_AM / p.R_p, rel=1e-15)


def test_condensed_step_solves_full_system():
    """At the condensed step's solution, a Newton update of the full (uncondensed) system is round-off,
    and the condensed iteration converges quadratically."""
    p = crystal(C_rate=1.0)
    m = CrystalModel(p)
    st = m.initial_state()
    for _ in range(20):
        st = m.newton_step(st, 1.0, p.i_app)
    hist = []
    new = m.newton_step(st, 1.0, p.i_app, history=hist)
    assert len(hist) <= 5 and hist[-1] < 1e-12
    J = m.jacobian_dense(new, st, 1.0, p.i_app)
    Re, Rc = m.residuals(new, st, 1.0, p.i_app)
    d = np.linalg.solve(J, -np.concatenate([Re.ravel(), Rc.ravel()]))
    assert np.abs(d).max() < 1e-12


def test_crystal_conservation():
    p = crystal(C_rate=1.0)
    m = CrystalModel(p)
    st = m.initial_state()
    x0 = st.copy()
    for _ in range(300):
        st = m.newton_step(st, 1.0, p.i_app)
    dx = m.mesh.dx[m.nodes]
    li = lambda s_: (m.cs_mean(s_) * p.vf_AM * dx).sum()
    q = p.i_app * 300 / p.F
    assert abs((li(st) - li(x0)) / q - 1.0) < 1e-11           # a discharge (I > 0) lithiates: + It/F
    salt = lambda s_: (m.el.eps_node * m.mesh.dx * m.conc(s_)).sum()
    assert abs(salt(st) / salt(x0) - 1.0) < 1e-12


def test_condensed_jacobian_matches_fd():
    """Full (uncondensed) Jacobian of (Re, Rc) vs central differences at a mid-discharge state."""
    p = crystal(C_rate=1.0, nj=26, sep_node=7, nj_crystal=6)
    m = CrystalModel(p)
    st = m.initial_state()
    for _ in range(200):
        st = m.newton_step(st, 1.0, p.i_app)
    old = st.copy()
    rng = np.random.default_rng(2)
    st = CrystalState(st.x + 1e-5 * rng.standard_normal(st.x.shape), st.xc + 1e-3 * rng.standard_normal(st.xc.shape))
    J = m.jacobian_dense(st, old, 1.0, p.i_app)
    flat = lambda s_: np.concatenate([a.ravel() for a in m.residuals(s_, old, 1.0, p.i_app)])
    z0 = np.concatenate([st.x.ravel(), st.xc.ravel()])
    ne = st.x.size
    Jfd = np.zeros_like(J)
    for q in range(z0.size):
        h = 1e-7 * max(1.0, abs(z0[q]))
        zp, zm = z0.copy(), z0.copy(); zp[q] += h; zm[q] -= h
        unpack = lambda z: CrystalState(z[:ne].reshape(st.x.shape), z[ne:].reshape(st.xc.shape))
        Jfd[:, q] = (flat(unpack(zp)) - flat(unpack(zm))) / (2 * h)
    scale = np.abs(Jfd).max(axis=1, keepdims=True) + 1e-300
    assert (np.abs(J - Jfd) / scale).max() < 1e-5


@pytest.mark.skip(reason="needs the protocol integration of Task 4")
@pytest.mark.parametrize("D_c", [1e-9, 1e-6])
def test_uniform_limit(D_c):
    """Fast solid diffusion: the crystal model's voltage approaches the uniform-particle model's."""
    from lfp_model.simulate import run
    a_u = run(Params(mode="corrected", C_rate=1.0), max_steps=1800).array
    a_c = run(crystal(C_rate=1.0, D_c=D_c), max_steps=1800).array
    n = min(len(a_u), len(a_c))
    err = np.abs(a_c[:n, 1] - a_u[:n, 1]).max()
    assert err < (2e-3 if D_c == 1e-9 else 5e-6), err
