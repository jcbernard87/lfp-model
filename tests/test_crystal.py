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
