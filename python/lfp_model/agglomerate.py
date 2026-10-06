"""Agglomerate scale: porous agglomerates of crystals, coupled to the electrode (docs/model.md section 13).

Corrected mode with particle_model = 'agglomerate'. A double-porosity model with one reaction description,
ported from nmc111-model's corrected agglomerate model with LiFePO4's kinetics:

- **Electrode scale** (macro-pores between agglomerates): the electrolyte concentration, phi1 and phi2. It has
  no reaction term of its own: the agglomerates at each cathode volume draw salt, ionic current and electronic
  current through their surfaces.
- **Agglomerate scale** (porous spheres of radius R_agg of uniform crystals of radius R_p): the pore
  electrolyte, phi1, phi2 and the crystal lithiation along the radius, with Butler-Volmer kinetics on the
  crystals, and the electrode's values at the surface.
- **Fully coupled Newton**: each Newton iteration condenses the agglomerates onto the electrode's diagonal
  blocks, with one factorization of all agglomerate systems stacked into one block-tridiagonal system.

The unknowns are the log variables of lfp_model.logcore. The electrode's fourth unknown has no equation of its
own here and is held fixed. A flux `q_in` is into an agglomerate through its surface.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional

import numpy as np

import bandsolver

from . import kinetics
from .driver import SolverFailure, limit_reason
from .logcore import (DIVERGED, N, P1, P2, S, U, Electrode, LogKinetics, Transport, bounded, converged, equilibrate,
                      logit, physical_update, sigmoid)
from .model import make_mesh
from .params import Params


@dataclass
class AggState:
    c: np.ndarray        # electrode (nj, 4): u, phi1, phi2, (fixed)
    ca: np.ndarray       # agglomerates at the cathode volumes (nl, nja, 4): u, phi1, phi2, s

    def copy(self):
        return AggState(self.c.copy(), self.ca.copy())


class AgglomerateModel:
    """Residuals, Jacobians and the condensed Newton step of the agglomerate model."""

    def __init__(self, p: Params):
        if p.particle_model != "agglomerate":
            raise ValueError("AgglomerateModel needs particle_model='agglomerate'")
        self.p = p
        self.mesh = make_mesh(p)
        m = self.mesh
        self.nodes = np.arange(m.s + 1, m.nj - 1)           # cathode volumes carry agglomerates
        self.nl = len(self.nodes)
        self.cs_max = kinetics.cs_max(p)
        self.kin = LogKinetics(coeffs=kinetics._ocp_coeffs(p), cs_max=self.cs_max, k=p.k_rxn, alpha_a=p.alpha_a,
                               alpha_c=p.alpha_c, c_bulk=p.c_bulk, R=p.R, T=p.T, F=p.F)
        self.tr = Transport(D0=p.D, t_plus=p.t_plus, c_bulk=p.c_bulk, kappa_bg=p.kappa_bg, R=p.R, T=p.T, F=p.F)
        self.el = Electrode(m, eps=p.eps, eps_sep=p.eps_sep, tau=p.tortuosity, tau_sep=p.tau_sep,
                            sigma=p.sigma, transport=self.tr)   # particles=None: the S column is frozen
        # radial grid (zero-volume center and surface nodes)
        na = p.nja
        h = p.R_agg / float(na - 2)
        j = np.arange(na)
        xa = np.zeros(na)
        xa[1:na - 1] = h * j[1:na - 1] - h / 2.0
        xa[na - 1] = p.R_agg
        dxa = np.zeros(na)
        dxa[1:na - 1] = h
        self.xa, self.dxa = xa, dxa
        r_face = xa[:-1] + dxa[:-1] / 2.0                   # face k is on the east side of node k
        r_face[0] = 0.0
        self.a_area = 4.0 * math.pi * r_face ** 2
        a_h = (dxa[:-1] + dxa[1:]) / 2.0
        r_W, r_E = xa - dxa / 2.0, xa + dxa / 2.0
        r_W[0] = r_E[0] = 0.0
        r_W[-1] = r_E[-1] = p.R_agg
        self.dV = 4.0 * math.pi / 3.0 * (r_E ** 3 - r_W ** 3)
        tau = p.tau_agg if p.tau_agg > 0 else p.eps_agg ** -0.5
        sig = p.sigma_agg if p.sigma_agg > 0 else p.sigma
        self.a_g = p.eps_agg / tau / a_h                     # per agglomerate face
        self.a_gs = (1.0 - p.eps_agg) * sig / a_h
        self.a_x = 3.0 * (1.0 - p.eps_agg) / p.R_p             # crystal area per agglomerate volume
        self.v_agg = p.vf_AM / (1.0 - p.eps_agg)                # agglomerate volume fraction
        self.s_agg = 3.0 * self.v_agg / p.R_agg                 # agglomerate surface per electrode volume

    # ------------------------------------------------------------------ state helpers
    def initial_state(self) -> AggState:
        p = self.p
        c = np.zeros((p.nj, N))
        c[:, P1], c[:, P2] = p.phi1_init, p.phi2_init
        c[:, S] = logit(p.cs_init / self.cs_max)
        ca = np.zeros((self.nl, p.nja, N))
        ca[..., P1], ca[..., P2] = p.phi1_init, p.phi2_init
        ca[..., S] = logit(p.cs_init / self.cs_max)
        return AggState(c, ca)

    def electrode(self, st: AggState) -> np.ndarray:
        return st.c

    def conc(self, st: AggState):
        """Electrolyte concentration [mol/cm3] of the electrode (nj,)."""
        return self.p.c_bulk * np.exp(st.c[:, U])

    def cs(self, st: AggState) -> np.ndarray:
        """Crystal concentrations [mol/cm3] (nl, nja)."""
        return self.cs_max * sigmoid(st.ca[..., S])

    # ------------------------------------------------------------------ agglomerate scale
    def agglomerates(self, xa, xaold, dt, xe):
        """Residual and blocks of every agglomerate (nl, na, ...); surface values xe (nl, 4)."""
        p = self.p
        nl, na = xa.shape[:2]
        ea, F = p.eps_agg, p.F
        Fv, dFa, dFb = self.tr.fluxes(xa[:, :-1], xa[:, 1:], self.a_g, self.a_gs)
        i, di = self.kin.rate(xa)
        Ar = self.a_area
        R = np.zeros_like(xa)
        A = np.zeros((nl, na, N, N)); B = np.zeros_like(A); D = np.zeros_like(A)
        for col in (U, P1, P2):                              # center: zero gradient
            R[:, 0, col] = xa[:, 1, col] - xa[:, 0, col]
            B[:, 0, col, col] = -1.0
            D[:, 0, col, col] = 1.0
        ji = np.arange(1, na - 1)                            # interior volumes
        V = self.dV[ji]
        for f, row in enumerate((U, P1, P2)):
            R[:, ji, row] = Ar[ji] * Fv[:, ji, f] - Ar[ji - 1] * Fv[:, ji - 1, f]
            B[:, ji, row] += Ar[ji, None] * dFa[:, ji, f] - Ar[ji - 1, None] * dFb[:, ji - 1, f]
            D[:, ji, row] += Ar[ji, None] * dFb[:, ji, f]
            A[:, ji, row] += -Ar[ji - 1, None] * dFa[:, ji - 1, f]
        c = p.c_bulk * np.exp(xa[..., U])
        cold = p.c_bulk * np.exp(xaold[..., U])
        R[:, ji, U] += ea * V * (c[:, ji] - cold[:, ji]) / dt - self.a_x * i[:, ji] * V / F
        B[:, ji, U, U] += ea * V * c[:, ji] / dt
        B[:, ji, U, :] -= self.a_x * di[:, ji] * V[None, :, None] / F
        R[:, ji, P1] += self.a_x * i[:, ji] * V
        B[:, ji, P1, :] += self.a_x * di[:, ji] * V[None, :, None]
        R[:, ji, P2] -= self.a_x * i[:, ji] * V
        B[:, ji, P2, :] -= self.a_x * di[:, ji] * V[None, :, None]
        for col in (U, P1, P2):                              # surface: the electrode's u, phi1, phi2
            R[:, -1, col] = xa[:, -1, col] - xe[:, col]
            B[:, -1, col, col] = 1.0
        # crystals, every node: (1 - eps_agg) cs_max dtheta/dt = -a_x i/F
        th, thold = sigmoid(xa[..., S]), sigmoid(xaold[..., S])
        R[..., S] = (1.0 - ea) * self.cs_max * (th - thold) / dt + self.a_x * i / F
        B[..., S, :] = self.a_x * di / F
        B[..., S, S] += (1.0 - ea) * self.cs_max * th * sigmoid(-xa[..., S]) / dt
        return R, A, B, D

    def surface_flux(self, xa):
        """q_in (nl, 3) = (N+, i1, i2) into each agglomerate through r = R, and dq/dx at its last two nodes
        (nl, 3, 2, 4)."""
        Fv, dFa, dFb = self.tr.fluxes(xa[:, -2], xa[:, -1], self.a_g[-1], self.a_gs[-1])
        return -Fv, -np.stack([dFa, dFb], axis=2)

    # ------------------------------------------------------------------ Newton
    def residuals(self, st: AggState, old: AggState, dt, I):
        """Full (uncondensed) residuals, for tests: (electrode (nj, 4), agglomerates (nl, na, 4))."""
        Re = self.el.residual_and_blocks(st.c, old.c, dt, I)[0]
        q, _ = self.surface_flux(st.ca)
        Re[self.nodes, :3] += (self.s_agg * self.mesh.dx[self.nodes])[:, None] * q
        Ra = self.agglomerates(st.ca, old.ca, dt, st.c[self.nodes])[0]
        return Re, Ra

    def newton_step(self, old: AggState, dt: float, I: float, *, backend: str = "fortran",
                    start: Optional[AggState] = None, history: Optional[list] = None) -> AggState:
        """One backward-Euler step from `old`, solved by the condensed Newton iteration."""
        p = self.p
        nl, na = self.nl, p.nja
        st = (start or old).copy()
        w = self.s_agg * self.mesh.dx[self.nodes]
        prev = math.inf
        for _ in range(p.newton_max_iter):
            # agglomerates: factor once (all stacked); solve for the update and the three surface responses
            Ra, Aa, Ba, Da = self.agglomerates(st.ca, old.ca, dt, st.c[self.nodes])
            rhs = np.zeros((4, nl, na, N))
            rhs[0] = -Ra
            for k, col in enumerate((U, P1, P2), start=1):
                rhs[k, :, na - 1, col] = 1.0
            Aa, Ba, Da, *r = equilibrate(Aa, Ba, Da, *rhs)
            try:
                fac = bandsolver.factor(Aa.reshape(-1, N, N), Ba.reshape(-1, N, N), Da.reshape(-1, N, N),
                                        backend=backend)
                sol = [fac.solve(rk.reshape(-1, N)).reshape(nl, na, N) for rk in r]
            except (bandsolver.NonFiniteError, bandsolver.SingularBlockError) as e:
                raise SolverFailure(f"agglomerate solve: {e}") from e
            dxa0 = sol[0]
            Z = np.stack(sol[1:], axis=-1)       # (nl, na, 4, 3): response to the electrode's u, phi1, phi2
            # electrode, with the condensed agglomerate sources: R_e + w q, linearized in the agglomerates
            Re, Ae, Be, De = self.el.residual_and_blocks(st.c, old.c, dt, I)
            q, Jq = self.surface_flux(st.ca)
            Re[self.nodes, :3] += w[:, None] * q
            Jf = Jq.reshape(nl, 3, 2 * N)
            Be[self.nodes, :3, :3] += w[:, None, None] * np.einsum("lrk,lkc->lrc", Jf,
                                                                    Z[:, na - 2:].reshape(nl, 2 * N, 3))
            G = -Re
            G[self.nodes, :3] -= w[:, None] * np.einsum("lrk,lk->lr", Jf, dxa0[:, na - 2:].reshape(nl, 2 * N))
            Ae, Be, De, G = equilibrate(Ae, Be, De, G)
            try:
                dce = bandsolver.solve(Ae, Be, De, G, backend=backend)
            except (bandsolver.NonFiniteError, bandsolver.SingularBlockError) as e:
                raise SolverFailure(f"electrode solve: {e}") from e
            dca = dxa0 + np.einsum("lnkc,lc->lnk", Z, dce[self.nodes, :3])
            lam = bounded((dce, dca))
            # divergence is judged on the damped electrode update: near theta = 0 or 1 a linearized log-odds
            # update of the agglomerates is legitimately huge, and through the condensation it also inflates
            # the undamped electrode update; `bounded` scales both down
            raw = lam * float(np.max(np.abs(dce)))
            st = AggState(st.c + lam * dce, st.ca + lam * dca)
            upd = max(physical_update(st.c, dce, frozen_s=True), physical_update(st.ca, dca))
            if history is not None:
                history.append(upd)
            if not math.isfinite(raw) or raw > DIVERGED:
                break
            if lam == 1.0 and converged(upd, prev, p.newton_tol):
                return st
            prev = upd
        raise SolverFailure("Newton did not converge")

    def limit_reason(self, st: AggState):
        c_min = min(float(self.conc(st).min()), float((self.p.c_bulk * np.exp(st.ca[..., U])).min()))
        th = sigmoid(st.ca[..., S])
        return limit_reason(c_min, self.p.c_bulk, float(th.min()), float(th.max()))
