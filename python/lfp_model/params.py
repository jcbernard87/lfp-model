"""Model parameters.

Defaults are the values of the original research code (see docs/parameters.md).
`Params.faithful()` returns the values that program actually stored: several
constants were written as single-precision Fortran literals and are rounded to
float32 before being used in double-precision arithmetic (deviation D-4).
"""
from __future__ import annotations

from dataclasses import dataclass, fields, replace

import numpy as np

MODES = ("faithful", "corrected")


def f32(x: float) -> float:
    """Round to the nearest single-precision value, returned as a Python float."""
    return float(np.float32(x))


@dataclass(frozen=True)
class Params:
    # --- geometry and mesh ---
    L_cath: float = 24.0e-4        # cathode thickness [cm]
    L_sep: float = 25.0e-4         # separator thickness [cm]
    nj: int = 101                  # total nodes
    sep_node: int = 22             # separator/cathode interface node (1-based, as in the original)
    eps: float = 0.5               # cathode porosity
    eps_AM: float = 0.8            # active-material volume fraction
    eps_sep: float = 0.39          # separator porosity
    tau_sep: float = 4.0           # separator tortuosity
    bruggeman: float = -0.5        # cathode tortuosity = eps**bruggeman
    # --- electrolyte ---
    D: float = 2.0e-6              # salt diffusivity [cm2/s]
    t_plus: float = 0.25           # cation transference number
    c_bulk: float = 1.0e-3         # initial electrolyte concentration [mol/cm3]
    z_plus: float = 1.0
    z_minus: float = -1.0
    # --- active material ---
    sigma: float = 3.0e-3          # electronic conductivity [S/cm]
    M: float = 125.759             # molar mass [g/mol]
    rho: float = 3.6               # density [g/cm3]
    Q_th: float = 0.170            # theoretical capacity [Ah/g]
    R_p: float = 200.0e-7          # particle radius [cm]
    k_rxn: float = 1.0e-8 * 10 ** 0.966   # rate constant
    alpha_a: float = 0.5
    alpha_c: float = 0.5
    # --- lithium counter electrode (output only) ---
    k_Li: float = 1.0e-6
    c_Li_ref: float = 1.0e-3
    # --- constants ---
    R: float = 8.314
    T: float = 298.0
    F: float = 96485.0
    # --- operation ---
    C_rate: float = 1.0            # [1/h], positive = discharge
    phi1_init: float = 3.6         # [V]
    phi2_init: float = 0.0         # [V]
    cs_init: float = 1.0e-5        # [mol/cm3]
    t_max: float = 36000.0         # [s]
    n_steps: int = 36000
    # --- numerics ---
    fd_step: float = 1.0e-6        # absolute step of the finite-difference reaction derivatives
    mode: str = "faithful"

    def __post_init__(self):
        if self.mode not in MODES:
            raise ValueError(f"mode must be one of {MODES}, got {self.mode!r}")

    @classmethod
    def faithful(cls, **overrides) -> "Params":
        """Parameters exactly as stored by the original program (float32-rounded literals).

        A `C_rate` override is also rounded to float32, because the original run generator
        wrote it into the source as a single-precision literal (for example `C_rate = 0.1`).
        """
        if "C_rate" in overrides:
            overrides["C_rate"] = f32(overrides["C_rate"])
        p = cls(
            R=f32(8.314), c_bulk=f32(0.001), Q_th=f32(0.170), M=f32(125.759), rho=f32(3.6),
            phi1_init=f32(3.6), eps_sep=f32(0.39), eps_AM=f32(0.8),
            k_rxn=1.0e-8 * f32(10.0 ** f32(0.966)), mode="faithful",
        )
        return replace(p, **overrides)

    def with_(self, **changes) -> "Params":
        return replace(self, **changes)

    def as_dict(self) -> dict:
        return {f.name: getattr(self, f.name) for f in fields(self)}

    # --- derived quantities (evaluated in the same order as the original) ---
    @property
    def mol_vol(self) -> float:
        return self.rho / self.M

    @property
    def spec_a(self) -> float:
        """Specific interfacial area a = 3 eps_AM / R_p [1/cm]."""
        return 3 * self.eps_AM / self.R_p

    @property
    def tortuosity(self) -> float:
        return self.eps ** self.bruggeman

    @property
    def i_specific(self) -> float:
        """Applied specific current [A/g]."""
        return self.Q_th * self.C_rate

    @property
    def i_app(self) -> float:
        """Applied current density [A/cm2]."""
        return self.i_specific * self.L_cath * self.eps_AM * self.rho

    @property
    def dt(self) -> float:
        return self.t_max / float(self.n_steps)
