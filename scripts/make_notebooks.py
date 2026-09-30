"""Generate the tutorial notebooks (committed without outputs; CI executes them)."""
from pathlib import Path

import nbformat as nbf

HERE = Path(__file__).resolve().parents[1] / "notebooks"


def nb(cells):
    n = nbf.v4.new_notebook()
    n.cells = [nbf.v4.new_markdown_cell(c[1]) if c[0] == "md" else nbf.v4.new_code_cell(c[1]) for c in cells]
    n.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3", "language": "python"}
    return n


quickstart = [
("md", """# 1 · Quickstart

This notebook runs the LiFePO₄ half-cell model from Python. The model is a porous LiFePO₄ cathode against a lithium-foil counter electrode, with a separator between them. The same model is also available as a self-contained Fortran program (`fortran/lfp.f90`) and C++ program (`cpp/lfp.cpp`), which read the same input file.

All parameters are the model's defaults (see `docs/parameters.md`), for illustration only. They are not fitted to any particular cell."""),
("code", """import numpy as np
import matplotlib.pyplot as plt

from lfp_model.params import Params
from lfp_model.simulate import run"""),
("md", """## A constant-current discharge

`mode="corrected"` is the recommended model: it fixes the defects of the original research code listed in `docs/deviations.md`. The run stops at the lower voltage cutoff `V_min` (2.5 V by default)."""),
("code", """p = Params(mode="corrected", C_rate=1.0)
r = run(p)
a = r.array          # columns: t [h], V, equivalents, Li eta [mV], Li i0 [mA/cm2], c at x=0, I [mA/cm2], step
print(r.exit_reason, f"after {r.steps} time steps")

capacity = a[:, 2] * p.F / (p.M * 3.6)      # electron equivalents -> mAh/g
plt.plot(capacity, a[:, 1])
plt.xlabel("capacity [mAh/g]"); plt.ylabel("cell voltage [V]"); plt.title("1C discharge")
plt.show()"""),
("md", """## Faithful and corrected modes

`Params.faithful()` reproduces the original program exactly. That includes ending the discharge on a numerical failure (NaN) rather than at a cutoff, and a different lithium-electrode overpotential (deviation D-12), which accounts for most of the voltage difference here."""),
("code", """rf = run(Params.faithful(C_rate=1.0))
af = rf.array
print("faithful run ended with:", rf.exit_reason)
ok = np.isfinite(af[:, 1])
plt.plot(af[ok, 2] * p.F / (p.M * 3.6), af[ok, 1], label="faithful (original code)")
plt.plot(capacity, a[:, 1], label="corrected")
plt.xlabel("capacity [mAh/g]"); plt.ylabel("cell voltage [V]"); plt.legend(); plt.show()"""),
("md", """## Running from an input file

The input file `input/default.nml` is shared by all three implementations:

```sh
python -m lfp_model input/default.nml        # Python
build/fortran/lfp_f input/default.nml        # Fortran (cmake --build build)
build/cpp/lfp_cpp   input/default.nml        # C++
```

Each writes `Time_Voltage.txt`. From Python you can also load the file and change values:"""),
("code", """from pathlib import Path
from lfp_model.namelist import load

root = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
p_file, extra = load(root / "input" / "default.nml")
print(p_file.mode, p_file.C_rate, extra)
r2 = run(p_file.with_(mode="corrected", C_rate=2.0))
print(r2.exit_reason, f"final capacity {r2.array[-1, 2] * p.F / (p.M * 3.6):.1f} mAh/g")"""),
]

explained = [
("md", r"""# 2 · The model, explained

This notebook walks through the equations (full detail in `docs/model.md`), the finite-volume grid, and one time step solved with **bandsolver**, Newman's BAND method.

## Unknowns and equations

At each node there are four unknowns: electrolyte concentration $c$, solid potential $\Phi_1$, electrolyte potential $\Phi_2$ and the lithium concentration in the particles $c_s$. In the porous cathode:

$$
\varepsilon\frac{\partial c}{\partial t} = -\frac{\partial N_+}{\partial x} + \frac{a\,i_n}{F},\qquad
\frac{\partial i_1}{\partial x} = -a\,i_n,\qquad
\frac{\partial i_2}{\partial x} = a\,i_n,\qquad
\varepsilon_{AM}\frac{\partial c_s}{\partial t} = -\frac{a\,i_n}{F}
$$

with the cation flux and the ionic current of a binary 1:1 electrolyte (dilute-solution theory)

$$
N_+ = -\varepsilon\frac{D_+}{\tau}\frac{\partial c}{\partial x} - \varepsilon\frac{u_+}{\tau}F c\frac{\partial \Phi_2}{\partial x},\qquad
i_2 = -\varepsilon F\frac{D_+ - D_-}{\tau}\frac{\partial c}{\partial x} - \varepsilon F^2\frac{u_+ + u_-}{\tau} c\frac{\partial\Phi_2}{\partial x},
$$

$i_1 = -(1-\varepsilon)\sigma\,\partial\Phi_1/\partial x$, and Butler–Volmer kinetics
$i_n = i_0\left[e^{\alpha_a F\eta/RT} - e^{-\alpha_c F\eta/RT}\right]$ with $\eta = \Phi_1 - \Phi_2 - U(\theta)$.
The separator has only transport. At the lithium foil $N_+ = I/F$ and $\Phi_2 = 0$; at the current collector $i_1 = I$.

## Log variables (corrected mode)

The corrected mode solves for $u = \ln(c/c_\text{bulk})$ and the particles' log-odds $s = \ln(\theta/(1-\theta))$, $\theta = c_s/c_{s,\max}$, instead of $c$ and $c_s$. Concentrations then stay positive and $0 < \theta < 1$ by construction, so the electrolyte can run out and particles can fill or empty smoothly, without clipping. The ion fluxes use exponential fitting (Scharfetter–Gummel), and the open-circuit potential has thermodynamic tails at the ends of the lithiation range (`docs/model.md` section 10). `LogModel.conc` and `LogModel.cs` convert a state back to concentrations."""),
("code", """import numpy as np
import matplotlib.pyplot as plt
import bandsolver

from lfp_model.logcore import logit, sigmoid
from lfp_model.logmodel import LogModel
from lfp_model.params import Params

p = Params(mode="corrected", C_rate=1.0)
model = LogModel(p)
m = model.mesh
print(f"{p.nj} nodes; separator/cathode interface at node {m.s}")"""),
("md", """## The grid

Control volumes of width $\\Delta x_j$, with zero-volume nodes at the lithium foil, at the separator/cathode interface and at the current collector. Those three nodes carry the boundary and interface conditions."""),
("code", """plt.figure(figsize=(8, 1.6))
plt.plot(m.x * 1e4, np.zeros_like(m.x), "|", ms=12)
plt.axvline(p.L_sep * 1e4, color="k", lw=0.8)
plt.xlabel("x [µm]"); plt.yticks([]); plt.title("node positions (separator | cathode)"); plt.show()"""),
("md", """## The open-circuit potential

The fitted LiFePO₄ curve is flat over most of the range and bounded at both ends. The thermodynamic tails (deviation D-16) add the ideal-solution terms $\\pm(RT/F)\\ln\\theta$ that dominate only very close to an empty or full particle (outside $10^{-4} < \\theta < 1 - 10^{-4}$), so a particle approaches either end asymptotically. Plotted against $s$, the tails are the straight sections at both ends."""),
("code", """s = np.linspace(-20, 20, 801)
U = model.kin.ocp(np.zeros_like(s), s)
fig, ax = plt.subplots(1, 2, figsize=(10, 3.5))
ax[0].plot(sigmoid(s), U); ax[0].set_xlabel("θ"); ax[0].set_ylabel("U [V]")
ax[1].plot(s, U); ax[1].set_xlabel("s = ln(θ/(1-θ))"); ax[1].set_ylabel("U [V]")
for a_ in ax: a_.set_ylim(2.8, 4.2)
plt.tight_layout(); plt.show()"""),
("md", """## One time step with BAND

Each backward-Euler step is a nonlinear system in block-tridiagonal form, $A_j\\,\\delta_{j-1} + B_j\\,\\delta_j + D_j\\,\\delta_{j+1} = G_j$ with $G = -R(x)$. The electrode assembles the residual and the blocks, and `bandsolver.solve` solves them. `newton_step` repeats this until the update is below the tolerance."""),
("code", """x0 = model.initial_state()
R, A, B, D = model.el.residual_and_blocks(x0, x0, 1.0, p.i_app)
print("block shapes:", A.shape, B.shape, D.shape, R.shape)
dx = bandsolver.solve(A, B, D, -R)
print("first-iteration update, max |dΦ1| =", np.abs(dx[:, 1]).max(), "V")

history = []
x1 = model.newton_step(x0, 1.0, p.i_app, history=history)
print("Newton converged in", len(history), "iterations; updates:", ", ".join(f"{h:.1e}" for h in history))"""),
("md", """## Checking the Jacobian

`bandsolver.check_jacobian` compares the hand-written Jacobian blocks with finite differences of the residual. The residual here includes the time term of the step from `x0`."""),
("code", """def fill(x):
    R, A, B, D = model.el.residual_and_blocks(x, x0, 1.0, p.i_app)
    return A, B, D, -R

chk = bandsolver.check_jacobian(fill, x1)
print(f"max relative error {chk.max_error:.1e}")"""),
("md", """## Inside the cell during a discharge"""),
("code", """from lfp_model.simulate import run

snapshots = {}
for t_s in (60, 1200, 2400, 3300):
    snapshots[t_s] = run(p, max_steps=t_s).final_state
fig, ax = plt.subplots(1, 3, figsize=(13, 3.5))
for t_s, x in snapshots.items():
    ax[0].plot(m.x * 1e4, model.conc(x) * 1e3, label=f"{t_s} s")
    ax[1].plot(m.x[m.s:] * 1e4, x[m.s:, 2] * 1e3)
    ax[2].plot(m.x[m.s:] * 1e4, model.cs(x)[m.s:] / model.cs_max)
ax[0].set_ylabel("c [M]"); ax[1].set_ylabel("Φ₂ [mV]"); ax[2].set_ylabel("θ")
for a_ in ax: a_.set_xlabel("x [µm]")
ax[0].legend(); plt.tight_layout(); plt.show()"""),
]

rates = [
("md", """# 3 · Rate capability and a full cycle

Discharges at several C-rates, and then a cycling protocol (docs/protocol.md) with constant-current and constant-voltage steps. These are illustrative default parameters, not fitted to a cell."""),
("code", """import numpy as np
import matplotlib.pyplot as plt

from lfp_model.params import Params
from lfp_model.simulate import run

to_mAhg = lambda p, eq: eq * p.F / (p.M * 3.6)"""),
("code", """rates = [0.2, 0.5, 1.0, 2.0, 5.0]
caps = []
for cr in rates:
    p = Params(mode="corrected", C_rate=cr)
    a = run(p).array
    plt.plot(to_mAhg(p, a[:, 2]), a[:, 1], label=f"{cr:g}C")
    caps.append(to_mAhg(p, a[-1, 2]))
plt.xlabel("capacity [mAh/g]"); plt.ylabel("cell voltage [V]"); plt.legend(); plt.show()

plt.semilogx(rates, caps, "o-")
plt.xlabel("C-rate"); plt.ylabel("capacity to 2.5 V [mAh/g]"); plt.show()"""),
("md", """## A full cycle

Discharge at 1C to 2.5 V, rest, charge at 1C to 4.0 V, hold 4.0 V until the current falls below C/20, then rest. The same protocol string works in the Fortran and C++ programs (`&protocol steps = '...'`)."""),
("code", """steps = "cc C=1 Vmin=2.5; rest t=1800; cc C=-1 Vmax=4.0; cv V=4.0 Imin=0.05; rest t=1800"
p = Params(mode="corrected", steps=steps)
r = run(p)
a = r.array
print(r.exit_reason)
fig, ax = plt.subplots(2, 1, figsize=(8, 5), sharex=True)
ax[0].plot(a[:, 0], a[:, 1]); ax[0].set_ylabel("V [V]")
ax[1].plot(a[:, 0], a[:, 6]); ax[1].set_ylabel("I [mA/cm²]"); ax[1].set_xlabel("time [h]")
plt.tight_layout(); plt.show()"""),
]

crystals = [
("md", """# 4 · The crystal scale

With `particle_model="crystal"` every cathode volume holds a crystal in which lithium diffuses (`docs/model.md` section 12), instead of a particle of uniform concentration. The crystal is a sphere by default (or a cylinder or slab), with radius `R_p` and solid diffusivity `D_c`. The crystals and the electrode are solved together in one Newton iteration per time step. Parameters are illustrative defaults, not fitted to a cell."""),
("code", """import numpy as np
import matplotlib.pyplot as plt

from lfp_model.crystal import CrystalModel
from lfp_model.params import Params
from lfp_model.simulate import run

to_mAhg = lambda p, eq: eq * p.F / (p.M * 3.6)
crystal = lambda **kw: Params(mode="corrected", particle_model="crystal", **kw)"""),
("md", """## Uniform particles and crystals

At the default D_c = 8 × 10⁻¹⁴ cm²/s, the diffusion time R_p²/D_c of a 200 nm crystal is about 5000 s, comparable to a 1C discharge, so solid diffusion costs capacity."""),
("code", """pu, pc = Params(mode="corrected", C_rate=1.0), crystal(C_rate=1.0)
au, ac = run(pu).array, run(pc).array
plt.plot(to_mAhg(pu, au[:, 2]), au[:, 1], label="uniform particles")
plt.plot(to_mAhg(pc, ac[:, 2]), ac[:, 1], label="crystals, D_c = 8e-14 cm²/s")
plt.xlabel("capacity [mAh/g]"); plt.ylabel("cell voltage [V]"); plt.title("1C discharge"); plt.legend(); plt.show()"""),
("md", """## Inside the crystals

Lithium profiles θ(r) after 1800 s at 1C, in crystals near the separator, in the middle of the electrode, and near the current collector. Each crystal fills from its surface (θ ≈ 0.6 at the surface, 0.36 at the centre). The three profiles nearly coincide: in this thin electrode, transport in the electrolyte and the solid phase is fast, and diffusion inside the crystals is what limits."""),
("code", """p = crystal(C_rate=1.0)
m = CrystalModel(p)
st = run(p, max_steps=1800).final_state
theta = m.cs(st) / m.cs_max
r_nm = m.cm.r * 1e7
for l, where in ((0, "near the separator"), (m.nl // 2, "middle"), (m.nl - 1, "near the collector")):
    plt.plot(r_nm, theta[l], "o-", ms=3, label=where)
plt.xlabel("r [nm]"); plt.ylabel("θ"); plt.title("crystal profiles after 1800 s at 1C"); plt.legend(); plt.show()"""),
("md", """## Capacity against rate and solid diffusivity"""),
("code", """rates = [1.0, 2.0, 5.0]
for D_c in (1e-12, 8e-14, 1e-14):
    caps = [to_mAhg(crystal(), run(crystal(C_rate=cr, D_c=D_c)).array[-1, 2]) for cr in rates]
    plt.semilogx(rates, caps, "o-", label=f"D_c = {D_c:g} cm²/s")
plt.xlabel("C-rate"); plt.ylabel("capacity to 2.5 V [mAh/g]"); plt.legend(); plt.show()"""),
("md", """## The uniform limit

As D_c grows the profiles flatten and the crystal model approaches the uniform-particle model: the voltage difference falls like 1/D_c."""),
("code", """a_u = run(Params(mode="corrected", C_rate=1.0), max_steps=900).array
Ds = [1e-12, 1e-11, 1e-10, 1e-9]
dV = []
for D_c in Ds:
    a_c = run(crystal(C_rate=1.0, D_c=D_c), max_steps=900).array
    n = min(len(a_u), len(a_c))
    dV.append(np.abs(a_c[:n, 1] - a_u[:n, 1]).max())
plt.loglog(Ds, dV, "o-", label="crystal − uniform")
plt.loglog(Ds, dV[0] * Ds[0] / np.array(Ds), "k--", lw=0.8, label="∝ 1/D_c")
plt.xlabel("D_c [cm²/s]"); plt.ylabel("max |ΔV| over 900 s at 1C [V]"); plt.legend(); plt.show()"""),
("md", """## Crystal shape

At the same R_p a sphere has three times, and a cylinder twice, the surface per volume of a slab, and a shorter mean diffusion path."""),
("code", """for shape in ("sphere", "cylinder", "slab"):
    p = crystal(C_rate=2.0, crystal_shape=shape)
    a = run(p).array
    plt.plot(to_mAhg(p, a[:, 2]), a[:, 1], label=shape)
plt.xlabel("capacity [mAh/g]"); plt.ylabel("cell voltage [V]"); plt.title("2C discharge"); plt.legend(); plt.show()"""),
]

for name, cells in (("01_quickstart", quickstart), ("02_model_explained", explained), ("03_rate_capability", rates),
                    ("04_crystal_scale", crystals)):
    nbf.write(nb(cells), HERE / f"{name}.ipynb")
print("wrote notebooks")
