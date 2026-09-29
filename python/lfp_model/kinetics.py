"""Open-circuit potential and Butler-Volmer kinetics (docs/model.md section 4).

All functions are vectorized over numpy arrays. In faithful mode the fit
coefficients and the literal 3.6 are float32-rounded, and the exchange current
density is rounded to float32 on every evaluation, as in the original program
(deviation D-4).
"""
from __future__ import annotations

import numpy as np

from .params import Params, f32

# Two-arctangent OCP fit: U = a0 + a1*atan(-b1*theta + c1) - a2*atan(-b2*theta + c2)
_OCP = (3.114559, 4.438792, 71.7352, 70.85337, 4.240252, 68.5605, 67.730082)


def _ocp_coeffs(p: Params):
    return tuple(f32(v) for v in _OCP) if p.mode == "faithful" else _OCP


def _lit36(p: Params) -> float:
    return f32(3.6) if p.mode == "faithful" else 3.6


def cs_max(p: Params) -> float:
    """Maximum lithium concentration in the active material [mol/cm3]."""
    return p.mol_vol * p.M * p.Q_th * 1000.0 * _lit36(p) / p.F


def theta(p: Params, cs):
    """State of lithiation cs/cs_max, evaluated as in the original (via x in Li_xFePO4)."""
    x = cs / p.mol_vol
    x_max = p.M * p.Q_th * 1000.0 * _lit36(p) / p.F
    return x / x_max


def ocp(p: Params, cs, c=None):
    """Open-circuit potential [V].

    Corrected mode adds the Nernst term (RT/F) ln(c/c_bulk) of the electrolyte concentration,
    which belongs in U with the electrostatic phi2 of the Nernst-Planck equations (deviation
    D-14). Faithful mode, and calls without c, return the original arctangent fit alone.
    """
    a0, a1, b1, c1, a2, b2, c2 = _ocp_coeffs(p)
    th = theta(p, cs)
    u = a0 + a1 * np.arctan(-(b1 * th) + c1) - a2 * np.arctan(-(b2 * th) + c2)
    if p.mode != "faithful" and c is not None:
        with np.errstate(divide="ignore", invalid="ignore"):
            u = u + p.R * p.T / p.F * np.log(c / p.c_bulk)
    return u


THETA_REG = 1.0e-6   # corrected mode: below this lithiation (or vacancy) fraction, x^alpha is regularized


def _power_reg(x, alpha, delta):
    """x**alpha for x >= delta; below it the C1 quadratic delta**alpha*((2-alpha)u + (alpha-1)u**2), u = x/delta.

    It keeps g(0) = 0 with a finite slope, so Newton's method stays well posed when a particle
    empties or fills completely (deviation D-13). Returns (g, dg/dx).
    """
    x = np.asarray(x, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        u = x / delta
        lo = x < delta
        g = np.where(lo, delta ** alpha * ((2.0 - alpha) * u + (alpha - 1.0) * u * u), x ** alpha)
        dg = np.where(lo, delta ** (alpha - 1.0) * ((2.0 - alpha) + 2.0 * (alpha - 1.0) * u), alpha * x ** (alpha - 1.0))
    return g, dg


def exchange_current(p: Params, c, cs):
    """i0 = F k c^aa (cs_max - cs)^aa cs^ac [A/cm2]."""
    if p.mode == "faithful":
        with np.errstate(invalid="ignore"):
            i0 = p.F * p.k_rxn * (c ** p.alpha_a) * ((cs_max(p) - cs) ** p.alpha_a) * (cs ** p.alpha_c)
        return np.asarray(i0, dtype=np.float32).astype(np.float64)
    return exchange_current_and_slope(p, c, cs)[0]


def exchange_current_and_slope(p: Params, c, cs):
    """Corrected-mode i0 and d(i0)/dcs, with the solid-concentration powers regularized (D-13)."""
    delta = THETA_REG * cs_max(p)
    gv, dgv = _power_reg(cs_max(p) - cs, p.alpha_a, delta)
    gs, dgs = _power_reg(cs, p.alpha_c, delta)
    with np.errstate(invalid="ignore"):
        pre = p.F * p.k_rxn * (c ** p.alpha_a)
    return pre * gv * gs, pre * (gs * -dgv + gv * dgs)


def reaction_rate(p: Params, c, cs, phi1, phi2):
    """Butler-Volmer current density per unit interfacial area [A/cm2]; anodic positive."""
    eta = phi1 - phi2 - ocp(p, cs, c)
    i0 = exchange_current(p, c, cs)
    rt = p.R * p.T
    with np.errstate(invalid="ignore", over="ignore"):
        return i0 * (np.exp(p.alpha_a * p.F * eta / rt) - np.exp(-(p.alpha_c * p.F * eta / rt)))


def reaction_derivatives(p: Params, c, cs, phi1, phi2):
    """Rate and its finite-difference derivatives (d/dc, d/dcs, d/dphi1, d/dphi2).

    Central differences with an absolute step, except forward differences when
    c or cs is at or below the step (as in the original, deviation D-6).
    """
    h = p.fd_step
    f = lambda c_, cs_, p1_, p2_: reaction_rate(p, c_, cs_, p1_, p2_)
    i = f(c, cs, phi1, phi2)
    with np.errstate(invalid="ignore"):
        d_c = np.where(c <= h,
                       (f(c + h, cs, phi1, phi2) - i) / h,
                       (f(c + h, cs, phi1, phi2) - f(c - h, cs, phi1, phi2)) / (2.0 * h))
        d_cs = np.where(cs <= h,
                        (f(c, cs + h, phi1, phi2) - i) / h,
                        (f(c, cs + h, phi1, phi2) - f(c, cs - h, phi1, phi2)) / (2.0 * h))
        d_p1 = (f(c, cs, phi1 + h, phi2) - f(c, cs, phi1 - h, phi2)) / (2.0 * h)
        d_p2 = (f(c, cs, phi1, phi2 + h) - f(c, cs, phi1, phi2 - h)) / (2.0 * h)
    return i, d_c, d_cs, d_p1, d_p2


def ocp_slope(p: Params, cs):
    """dU/dcs [V cm3/mol]."""
    a0, a1, b1, c1, a2, b2, c2 = _ocp_coeffs(p)
    th = theta(p, cs)
    dth = 1.0 / (p.mol_vol * (p.M * p.Q_th * 1000.0 * _lit36(p) / p.F))
    du = a1 * (-b1) / (1.0 + (-(b1 * th) + c1) ** 2) - a2 * (-b2) / (1.0 + (-(b2 * th) + c2) ** 2)
    return du * dth


def reaction_derivatives_analytic(p: Params, c, cs, phi1, phi2):
    """Rate and exact derivatives (d/dc, d/dcs, d/dphi1, d/dphi2); used in corrected mode (fixes D-6)."""
    rt = p.R * p.T
    A_, B_ = p.alpha_a * p.F / rt, p.alpha_c * p.F / rt
    eta = phi1 - phi2 - ocp(p, cs, c)
    i0, di0_dcs = exchange_current_and_slope(p, c, cs)
    with np.errstate(invalid="ignore", over="ignore", divide="ignore"):
        ea, ec = np.exp(A_ * eta), np.exp(-B_ * eta)
        i = i0 * (ea - ec)
        di_deta = i0 * (A_ * ea + B_ * ec)
        d_c = p.alpha_a * i / c - di_deta * (rt / p.F) / c        # includes dU/dc of the Nernst term
        d_cs = di0_dcs * (ea - ec) - di_deta * ocp_slope(p, cs)
    return i, d_c, d_cs, di_deta, -di_deta


def li_exchange_current(p: Params, c):
    """Exchange current density of the lithium counter electrode [A/cm2]."""
    return p.F * p.k_Li * (c ** 0.5) * (p.c_Li_ref ** 0.5)


def li_foil(p: Params, c, I):
    """The lithium foil's Nernst potential U_Li and overpotential eta_Li [V] (corrected mode).

    With the foil metal as the 0 V reference, phi2 at the foil is -(U_Li + eta_Li), where
    U_Li = (RT/F) ln(c/c_Li_ref) and eta_Li = (RT/(alpha F)) asinh(I/(2 i0)) (symmetric
    Butler-Volmer, alpha = 0.5; I > 0 on discharge, when the foil is oxidized).
    """
    alpha = 0.5
    rtf = p.R * p.T / p.F
    i0 = li_exchange_current(p, c)
    return rtf * np.log(c / p.c_Li_ref), rtf / alpha * np.arcsinh(I / (2.0 * i0))
