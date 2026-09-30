# LFP model: equations and discretization

This document describes the model **as implemented in the original research code**, which is the reference for the faithful port (`mode = "faithful"`). Suspected defects are marked **[D-n]** and tracked in [deviations.md](deviations.md). Source line numbers (`L123`) refer to the private reference copy of the original program (`LFP_Legacy.f95`, SHA-256 in the private oracle's `PROVENANCE.sha256`); the source itself is not distributed.

## 1. Cell and domain

A one-dimensional half cell: lithium-metal counter electrode | separator | porous LiFePO₄ cathode | current collector.

```
x = 0                x = L_sep                           x = L_sep + L_cath
 |---- separator ----|----------- porous cathode ----------|
node 1            SEP_NODE = 22                         NJ = 101
(Li foil face)    (separator/cathode interface)       (current collector)
```

- **Nodes (L890–L931).** A finite-volume grid with zero-volume boundary nodes at `j = 1`, `j = SEP_NODE` and `j = NJ`. Separator control volumes are `j = 2 … SEP_NODE−1` (20 volumes of width `h_sep = L_sep/(SEP_NODE−2)`). Cathode control volumes are `j = SEP_NODE+1 … NJ−1` (78 volumes of width `h_cath = L_cath/(NJ−SEP_NODE−1)`). Nodes sit at control-volume centres; `Δx_j = 0` at the three boundary nodes.
- **Face interpolation (L1119–L1491).** For the face between nodes `j` and `j+1`: `α_E = Δx_j/(Δx_j+Δx_{j+1})`, `β_E = 2/(Δx_j+Δx_{j+1})`, face value `c_E = α_E c_{j+1} + (1−α_E) c_j`, face gradient `(∂c/∂x)_E = β_E (c_{j+1} − c_j)`. The west face is analogous. With zero-volume boundary nodes this gives exact face values at the boundaries and central differences inside.

## 2. Unknowns (electrode scale, N = 4)

| k | Symbol | Meaning | Unit |
|---|---|---|---|
| 1 | c | electrolyte (salt) concentration, binary 1:1 electrolyte, c₊ = c₋ = c | mol/cm³ |
| 2 | Φ₁ | solid-phase (electronic) potential | V |
| 3 | Φ₂ | electrolyte potential | V |
| 4 | c_s | lithium concentration in the active material, per volume of active material | mol/cm³ |

A crystal (particle) scale with two unknowns exists in the source but is **inactive**: see §8.

## 3. Properties (L1029–L1050)

- Ion diffusivities from the salt diffusivity `D` and cation transference number `t₊` (anion `t₋ = 1 − t₊`), with `r = t₋/t₊`: `D₊ = D (1 + r)/(2r)`, `D₋ = r D₊`. These reproduce the binary salt diffusivity `2D₊D₋/(D₊+D₋) = D` and `t₊ = D₊/(D₊+D₋)`.
- Mobilities (Nernst–Einstein): `u_i = D_i/(RT)`.
- Effective properties: divide by tortuosity τ and multiply by porosity ε in the flux expressions. Cathode: `τ = ε^(−1/2)` (Bruggeman, so ε·D/τ = ε^1.5 D). Separator: `τ_sep = 4` with `ε_sep`.
- Electrolyte conductivity is not used as a separate parameter; the current is written directly from the ion fluxes (below).

**Fluxes** (per unit cross-section, ε and τ of the local region):

- Cation flux: `N₊ = −ε (D₊/τ) ∂c/∂x − ε z₊ (u₊/τ) F c ∂Φ₂/∂x`
- Ionic current: `i₂ = −ε F (z₊D₊ + z₋D₋)/τ ∂c/∂x − ε F² (z₊²u₊ + z₋²u₋)/τ c ∂Φ₂/∂x`
- Electronic current: `i₁ = −(1−ε) σ ∂Φ₁/∂x`

## 4. Kinetics (MOD_echem_rxn, L143–L271)

- **Butler–Volmer (L255–L269):** `i_n = i₀ [exp(α_a F η/RT) − exp(−α_c F η/RT)]`, with `η = Φ₁ − Φ₂ − U(θ)` and `α_a = α_c = 0.5`. Anodic (delithiation) current is positive.
- **Exchange current density (L227–L237):** `i₀ = F k c^α_a (c_s,max − c_s)^α_a c_s^α_c`.
- **Maximum solid concentration:** `c_s,max = ρ Q_th · 3600 / F`, from density and theoretical capacity (L233). The state of lithiation is `θ = c_s/c_s,max` (L156–L160).
- **Open-circuit potential (L151–L208):** the active expression is a two-arctangent fit: `U(θ) = 3.114559 + 4.438792·atan(−71.7352 θ + 70.85337) − 4.240252·atan(−68.5605 θ + 67.730082)`. A 31-term Redlich–Kister sum (`Vint`) is computed on every call and then discarded, and the concentration argument is unused. **[D-5]**
- **Reaction derivatives (L1060–L1083):** `∂i_n/∂c, ∂i_n/∂c_s, ∂i_n/∂Φ₁, ∂i_n/∂Φ₂` are found by finite differences with an absolute step of 10⁻⁶. They are central differences, except forward differences when `c ≤ 10⁻⁶` or `c_s ≤ 10⁻⁶`. A step of 10⁻⁶ mol/cm³ is about 4 % of `c_s,max` and 100 % of `c_s` at the start of discharge, so these derivatives are inaccurate. **[D-6]**

## 5. Governing equations

`a = 3 ε_AM / R_p` is the specific interfacial area (L127), with `R_p = xmax_c` (particle radius) and `ε_AM` the active-material volume fraction.

**Separator interior** (`2 ≤ j < SEP_NODE`, L1494–L1580):

| Row | Equation |
|---|---|
| 1 | `ε_sep ∂c/∂t = −∂N₊/∂x` |
| 2 | `∂/∂x[(1−ε_sep) σ ∂Φ₁/∂x] = 0` (a solid potential is carried in the separator but has no physical meaning there) |
| 3 | `∂i₂/∂x = 0` |
| 4 | `∂c_s/∂t = 0` |

**Cathode interior** (`SEP_NODE < j < NJ`, L1584–L1666):

| Row | Equation |
|---|---|
| 1 | `ε ∂c/∂t = −∂N₊/∂x + a i_n/F` |
| 2 | `∂i₁/∂x = −a i_n` |
| 3 | `∂i₂/∂x = +a i_n` |
| 4 | `ε_AM ∂c_s/∂t = −a i_n/F` (uniform particle: no solid-state diffusion) |

## 6. Boundary and interface conditions

**x = 0, Li-foil face (j = 1, L1119–L1207):**

| Row | Condition |
|---|---|
| 1 | cation flux into the cell equals the applied current: `N₊ = I/F`. **[D-2]** The flux is built with the cathode porosity `ε` and separator coefficients, not `ε_sep`. |
| 2 | intended zero electronic current, `∂Φ₁/∂x = 0`. **[D-1]** The sign of the correction term makes the linearized row double the old gradient instead of zeroing it (the same pattern as numerical-audit §8.19). |
| 3 | `Φ₂ = 0` (reference potential). The source notes that a current boundary condition "was not working". |
| 4 | `c_s` frozen (`∂c_s/∂t = 0`). |

**Separator/cathode interface (j = SEP_NODE, zero volume, L1218–L1357):**

| Row | Condition |
|---|---|
| 1 | cation-flux continuity; separator side built with `ε` (not `ε_sep`) and separator transport coefficients **[D-2]**; cathode side with `ε` and cathode coefficients |
| 2 | electronic-current continuity, `(1−ε)σ` on both sides |
| 3 | ionic-current continuity, with the same `ε` issue as row 1 **[D-2]** |
| 4 | `ε_AM ∂c_s/∂t = −a i_n/F` evaluated at the interface node (this node's `c_s` is used for the output OCP) |

**x = L, current collector (j = NJ, L1367–L1447):**

| Row | Condition |
|---|---|
| 1 | no cation flux, `N₊ = 0` |
| 2 | `i₁ = I` (applied current density, positive on discharge) |
| 3 | no ionic current, `i₂ = 0` |
| 4 | `ε_AM ∂c_s/∂t = −a i_n/F` at the boundary node |

## 7. Linearization and time stepping

- **Row form.** Every row is assembled as flux differences plus a local term: `F_W − F_E + Q_j(c) = 0`, where `F_W`, `F_E` are the face fluxes above and `Q_j` holds storage and reaction terms multiplied by `Δx_j` (rows 1–3). Row 4 is per unit volume and is not multiplied by `Δx_j`. Face fluxes are linearized about the previous state; `Q_j` is linearized with the finite-difference reaction derivatives. The block-tridiagonal system `A_j δ_{j−1} + B_j δ_j + D_j δ_{j+1} = G_j` (G = −residual) is solved once with BAND, and the state is updated `c ← c + δ` (L788–L824).
- **One linear solve per time step.** The time scheme is a linearized implicit (backward) Euler with no Newton iteration: the source comment (L769–L777) assumes the step is small enough for the linearization to be exact. **[D-7]**
- **Step size (L754–L760, L885).** `Δt = t_max/N_steps` (= 1 s with the stored `t_max = 36 000 s`, `N_steps = 36 000`), or `Δt ← 1.0001 Δt` during rest (`state = 'R'`, unused in these runs).
- **Coulomb counting (L996–L1009).** `mAh/g` is advanced by `1000·I_spec·Δt/3600` inside the `j = 1` fill call, before the solve.
- **Loop and exits (L685–L724).** Up to `N_steps` iterations. At the top of each iteration the program writes an output row (see §9) or exits:
  - end of charge: `Φ₁(NJ) ≥ 99` with `state = 'C'` (a sentinel, not a physical criterion)
  - `δc(1,1)` is NaN
  - `t ≥ 99 h`

  **There is no voltage cutoff:** every archived discharge ends on the NaN exit. **[D-3]**

## 8. Inactive code (present in the source, not executed in any archived run)

- **Crystal scale** (`crystal_scale = 0`, L104; routines `bound_val_c`, `fillmat_c`, and a duplicate `_c` solver). It is designed as solid diffusion in a slab of thickness `xmax_c` with an α/β two-phase model, coupled at the surface through `i_n`. It is not ported in faithful mode. Known problems if it is enabled:
  - `step_c_alpha` is undeclared and never set; `step_c_a` is the step actually used (L1716, L1726–L1755). **[D-8]**
  - `Rxn_Beta` tests an undeclared `theta_beta` instead of its argument `thet_beta` (L305). **[D-8]**
  - Two inconsistent saturation values: `c_alpha_sat = 1000` in `Rxn_Beta` vs `140000·ρ/M` in `fillmat_c`. **[D-8]**
  - The surface-flux coupling carries the comment "not sure this is right" (L1823–L1826).
- **Anode OCP** (`OCP_ANODE`) is unused. Its return value would be single precision.
- **Contact resistance** `R_contact` is never set; it is used only for an unused `eta_contact`.

## 9. Output (`write_all_voltage`, L380–L666)

`Time_Voltage.txt`: two header rows, then one row per write event with columns:

| Column | Header | Meaning |
|---|---|---|
| 1 | State | `D`/`C`/`R` |
| 2 | Time (hours) | `t/3600` |
| 3 | Voltage (Volts) | `Φ₁(NJ) + η_Li`, cathode potential relative to Φ₂(0) plus the lithium-anode overpotential |
| 4 | Equivalence | `mAh/g · M/(3.6 F)`, electrons per formula unit passed |
| 5 | Anode_Eta (mV) | `η_Li = ∓(RT/(2αF))·ln(I/i₀,Li)` on discharge/charge, with `α = 0.5` and `i₀,Li = F·10⁻⁶·c(0)^0.5·(10⁻³)^0.5` |
| 6 | anode_exchange_c (mA/cm²) | `i₀,Li` |
| 7 | Edge_c0 (mol/cm³) | `c` at x = 0 |

- The voltage excludes the concentration (Nernst) term of the lithium electrode. It is computed as `anode_potential` but only printed to the terminal. Corrected mode reports the voltage against the lithium foil, including this term **[D-15]**.
- **Write schedule.** One row at `it = 1` (t = 0, before the first solve); then whenever `(t − t_last)/3600 ≥ t_max/N_steps/200` (= 18 s); then one final row at the exit.
  - `t_last = t − Δt` is stored in an implicitly **integer** variable, so it is truncated. **[D-4]** With Δt = 1 s, rows are written every 17 s.
  - Rows are written *before* the step, so each row reports the state at time `t` with `mAh/g` already incremented to `t`.
- Format: `A5, 2F12.5, ES15.5…`.

## 10. Corrected mode: log variables, exponential fitting and OCP tails (v0.3.0)

Sections 1–9 describe the original equations in the concentrations c and c_s. Since v0.3.0 corrected mode solves the same conservation laws in variables that keep every state physical, so the electrolyte can run out and particles can fill or empty smoothly, without clipping. Faithful mode is unchanged.

**Unknowns.** At every node: u = ln(c/c_bulk), Φ₁, Φ₂ and the particles' log-odds s = ln(θ/(1−θ)), θ = c_s/c_s,max. Then c = c_bulk·eᵘ > 0 and 0 < θ < 1 by construction. The storage terms stay in the conserved quantities, ε(c − c_old)/Δt and ε_AM c_s,max(θ − θ_old)/Δt, so salt and lithium are conserved to round-off. θ is evaluated with a sigmoid that is accurate for either sign of s, and 1 − θ as sigmoid(−s), so neither rounds to zero.

**Kinetics.** The overpotential and the exchange current are written in the new variables:

  η = Φ₁ − Φ₂ − U_fit(θ) − U_tail(s) − (RT/F)·u
  ln i₀ = ln(F k c_bulk^α_a c_s,max^(α_a+α_c)) + α_a u − α_a ln(1 + eˢ) − α_c ln(1 + e⁻ˢ)

The electrolyte Nernst term (D-14) is linear in u. i₀ has bounded derivatives in (u, s) as θ → 0 or 1, which replaces the regularization D-13.

**Thermodynamic tails (D-16).** The arctangent fit U_fit(θ) is bounded at θ = 0 and θ = 1. With i₀ ∝ θ^α_c(1−θ)^α_a vanishing only as a power, a particle would then empty or fill completely in finite time, which no real insertion electrode does. The OCP gets the ideal-solution terms that dominate at the ends:

  U_tail(s) = (RT/F)·[ln(1 + e^(s_e − s)) − ln(1 + e^(s + s_e))],  s_e = ln(θ_e/(1 − θ_e)),  θ_e = 10⁻⁴

For θ ≪ θ_e this is (RT/F)·ln(θ_e/θ), and for 1 − θ ≪ θ_e it is −(RT/F)·ln(θ_e/(1−θ)): the potential diverges logarithmically and each end is approached asymptotically. Between 1 % and 99 % lithiation the tail changes U by at most 0.25 mV.

**Ion fluxes: exponential fitting (Scharfetter–Gummel).** For ion i with charge zᵢ across a face between nodes a and b, a distance h apart, with Δ = zᵢF(Φ₂,b − Φ₂,a)/(RT):

  Nᵢ = (ε/τ)(Dᵢ/h)·[B(Δ)·c_a − B(−Δ)·c_b],  B(x) = x/(eˣ − 1)

This is exact for a constant flux in a constant field between the nodes, keeps concentrations positive, and stays accurate where migration dominates (the depleted regime). For |Δ| ≪ 1 it reduces to the centred scheme of section 5. D₊ and D₋ are those of section 3. The cation row uses N₊, and the current row uses i₂ = F(N₊ − N₋), so charge and salt are consistent by construction.

**Background conductivity.** Where the salt is exhausted, no current can flow and Φ₂ is undefined. The solvent's own ionic conductivity κ_bg (`kappa_bg`, default 10⁻⁸ S/cm, about 10⁻⁶ of the 1 M electrolyte's) is added as an ohmic current −(ε/τ)κ_bg ∂Φ₂/∂x that carries no salt. It keeps Φ₂ defined there and changes nothing measurable elsewhere.

**Newton.**
- Each iteration's step is limited to |Δu| ≤ 1, |ΔΦ| ≤ 0.1 V and |Δs| ≤ 2. The exponentials make larger steps overshoot.
- Convergence is judged on the physical variables: max(c/c_bulk·|Δu|, |ΔΦ|, θ(1−θ)|Δs|) ≤ 10⁻¹⁰, or stagnation at the round-off floor (below 10⁻⁷ and down by less than half since the previous iteration). Near a limit the log variables are ill-conditioned, and their round-off is physically irrelevant.
- A raw update above 10³ means the step has no solution, for example a current that the remaining capacity cannot carry for the whole step. The driver then shortens the step (docs/protocol.md).

## 11. Branch checklist (T2 verification)

Every branch of `fillmat` (L980–L1670) is covered above:

- [x] Coulomb counting, `j = 1` (L996–L1009) → §7
- [x] property block (L1011–L1050) → §3
- [x] reaction rate and FD derivatives (L1057–L1083) → §4
- [x] `j = 1` rows 1–4 (L1119–L1207) → §6
- [x] `j = SEP_NODE` rows 1–4 (L1218–L1357) → §6
- [x] `j = NJ` rows 1–4 (L1367–L1447) → §6
- [x] separator interior rows 1–4 (L1494–L1580) → §5
- [x] cathode interior rows 1–4 (L1584–L1666) → §5
- [x] `fillmat_c` (L1679–L1939), not executed → §8

## 12. Crystal scale: solid diffusion (corrected mode, v0.4.0)

With `particle_model = 'crystal'` each cathode volume (electrode nodes between the interface and the collector) carries one representative crystal with solid-state diffusion, instead of a uniform particle. The original's crystal-scale code (§8) is not ported; this scale is written fresh. It is stage 1 of two: an α/β phase-change model of LiFePO₄ may follow in a later version on the same structure.

**Geometry.** `crystal_shape` = `slab`, `cylinder` or `sphere` (k = 0, 1, 2), with size R_p (the half-thickness of a slab, the radius otherwise). The crystal surface per electrode volume is a = (k+1)·ε_AM/R_p, which is the uniform model's 3ε_AM/R_p for a sphere.

**Equations.** With θ = c_s/c_s,max and a constant solid diffusivity D_c (`D_c`):

  c_s,max ∂θ/∂t = (1/rᵏ) ∂/∂r (rᵏ D_c c_s,max ∂θ/∂r),  ∂θ/∂r = 0 at r = 0,  −D_c c_s,max ∂θ/∂r = i_n/F at r = R_p (the outward flux; an anodic i_n removes lithium)

i_n is the Butler–Volmer rate of §10, evaluated with the electrode node's u, Φ₁, Φ₂ and the crystal's **surface** θ (OCP fit, tails D-16 and Nernst term D-14). The electrode rows carry the source a·i_n as in the uniform model; the electrode's own fourth unknown is held fixed.

**Discretization.** Vertex-centred finite volumes: nodes at r_j = j·h, h = R_p/(nj_crystal − 1), each owning the volume ∫ rᵏ dr over [r_j − h/2, r_j + h/2] ∩ [0, R_p]. The flux through the face between nodes j and j+1 is −D_c c_s,max A_f (θ_{j+1} − θ_j)/h with A_f = r_fᵏ; the surface node receives the reaction flux R_pᵏ·i_n/F. Every node, the surface included, stores lithium: a surface that is filling or emptying can take up or give off lithium, as in the continuum (a zero-volume surface node could pass only what diffusion across half a cell carries, which fails at small D_c or in a constant-voltage hold). The unknown at every crystal node is the log-odds s = ln(θ/(1−θ)), so 0 < θ < 1 throughout (§10). Lithium is conserved exactly: the crystal rows telescope to d/dt ∑ c_s,max V_j θ_j = −R_pᵏ i_n/F.

**Condensed Newton step.** The crystals couple only to their own electrode node, through the surface reaction. Every Newton iteration:
1. assembles all crystals as one block-tridiagonal system with blocks of size 1 (J_cc) and solves it for the residual and for three right-hand sides, the surface row's derivatives with respect to the node's u, Φ₁ and Φ₂ (J_ce);
2. eliminates the crystals: δc = −J_cc⁻¹(R_c + J_ce δe), whose surface entry adds a 3 × 3 term to each cathode node's diagonal block and a correction to its residual;
3. solves the electrode system with BAND and back-substitutes for δc.

This is the exact Newton step of the coupled system (quadratic convergence). The step limits of §10 apply to the crystal unknowns too. Divergence (an update above 10³) is judged on the electrode unknowns only: near θ = 0 or 1 a linearized log-odds update of the crystals is legitimately large, and the step limit damps it.

**Limits.** The exit reasons `particles_full` and `particles_empty` refer to the crystal surfaces.

**Uniform limit.** For D_c → ∞ the crystal profiles flatten and the model becomes the uniform-particle model; the difference in cell voltage falls like 1/D_c ([validation.md](validation.md) §6).
