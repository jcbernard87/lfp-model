"""Crystal scale: solid diffusion in the crystals, coupled to the electrode (docs/model.md section 12).

Corrected mode with particle_model = 'crystal'. Each interior cathode node carries one crystal
(slab, cylinder or sphere, k = 0, 1, 2), discretized by finite volumes in r with zero-volume
center and surface nodes. The unknown at every crystal node is the log-odds s = ln(theta/(1-theta)),
and the flux is Fick's law in theta with a constant D_c. Each Newton iteration condenses the
crystals into the electrode's diagonal blocks (the nmc111-model agglomerate pattern).
"""
from __future__ import annotations

import numpy as np

from .logcore import sigmoid


class CrystalMesh:
    def __init__(self, R: float, nc: int, k: int):
        h = R / float(nc - 2)
        j = np.arange(nc)
        r = np.zeros(nc)
        r[1:nc - 1] = h * j[1:nc - 1] - h / 2.0
        r[nc - 1] = R
        rW, rE = np.zeros(nc), np.zeros(nc)
        rW[1:nc - 1] = h * (j[1:nc - 1] - 1)
        rE[1:nc - 1] = h * j[1:nc - 1]
        self.r, self.k, self.R = r, k, R
        self.V = (rE ** (k + 1) - rW ** (k + 1)) / (k + 1)
        rf = h * np.arange(nc - 1)
        rf[-1] = R
        self.A = rf ** k                      # k = 0: 1 everywhere (face 0 is never used: symmetry)
        self.d = np.full(nc - 1, h)
        self.d[0] = self.d[-1] = h / 2.0
        self.Vtot = R ** (k + 1) / (k + 1)


def diffusion_rows(cm: CrystalMesh, D: float, cs_max: float, xc, xcold, dt):
    """Residual and tridiagonal dR/ds (sub An, diagonal B, super Dn) of the crystal rows, without the
    reaction. Every array is (nl, nc). The surface row is the outward flux per area; the caller
    subtracts i/F."""
    th, thold = sigmoid(xc), sigmoid(xcold)
    dth = th * sigmoid(-xc)
    g = D * cs_max * cm.A / cm.d                           # (nc-1,)
    Jf = -g * (th[:, 1:] - th[:, :-1])                     # outward flux through each face
    R = np.zeros_like(xc); B = np.zeros_like(xc); An = np.zeros_like(xc); Dn = np.zeros_like(xc)
    R[:, 0] = xc[:, 0] - xc[:, 1]
    B[:, 0] = 1.0; Dn[:, 0] = -1.0
    ji = np.arange(1, xc.shape[1] - 1)
    V = cm.V[ji]
    R[:, ji] = cs_max * V * (th[:, ji] - thold[:, ji]) / dt + Jf[:, ji]
    R[:, ji[1:]] -= Jf[:, ji[1:] - 1]                       # J_0 = 0: no flux through the center face
    B[:, ji] = cs_max * V * dth[:, ji] / dt + g[ji] * dth[:, ji]
    Dn[:, ji] = -g[ji] * dth[:, ji + 1]
    B[:, ji[1:]] += g[ji[1:] - 1] * dth[:, ji[1:]]
    An[:, ji[1:]] = -g[ji[1:] - 1] * dth[:, ji[1:] - 1]
    n = xc.shape[1] - 1
    R[:, n] = Jf[:, n - 1] / cm.A[n - 1]
    B[:, n] = -D * cs_max * dth[:, n] / cm.d[n - 1]
    An[:, n] = D * cs_max * dth[:, n - 1] / cm.d[n - 1]
    return R, B, Dn, An
