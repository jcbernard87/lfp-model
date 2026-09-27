# Deviations from the original research code

The faithful mode (`mode = "faithful"`) reproduces the original program, including every item below. The corrected mode (`mode = "corrected"`, the default once T8 is done) applies the fixes marked **fixed**. Each fix is a separate commit that states the effect on results.

Status values: **candidate** (suspected from reading the source), **confirmed** (demonstrated by a test or run), **fixed**, **kept** (reviewed and left as is, with the reason).

| ID | Status | Summary |
|---|---|---|
| D-1 | confirmed (latent) | Li-foil face (j = 1): the solid-potential row doubles the old gradient instead of zeroing it |
| D-2 | confirmed | Separator faces next to the boundary nodes use the cathode porosity `ε`: the separator transports as if ε = 0.5 while storing with ε_sep = 0.39, and 0.52 % of the salt is lost during start-up |
| D-3 | confirmed | No voltage cutoff: every discharge ends on a NaN |
| D-4 | confirmed | Precision and typing: single-precision literals, implicitly single-precision `ex_1`, implicitly integer `last_write_time` |
| D-5 | candidate | The OCP routine computes and discards a 31-term Redlich–Kister sum on every call |
| D-6 | candidate | Reaction derivatives use absolute FD steps of 10⁻⁶ (≈ 4 % of c_s,max, 100 % of the initial c_s) |
| D-7 | candidate | One linearized solve per time step (no Newton iteration); nonlinearity error is not controlled |
| D-8 | candidate | The inactive crystal-scale code has undeclared and misspelled variables and inconsistent constants |
| D-9 | candidate | Porosity 0.5 plus active-material fraction 0.8 add up to more than 1 |
| D-11 | confirmed | Interior and interface ionic-current rows leave the diffusion current out of the residual, so the model solves Ohm's law for Φ₂ and drops the diffusion potential |
| D-10 | confirmed | The output routine reads `cprev` at index (SEP_NODE−NJ)/2 = −39 (out of bounds; value unused) |

## D-1. Li-foil face solid-potential row sign

- **Where:** `fillmat`, j = 1, row 2 (L1168–L1169).
- **What:** the row is intended to impose `i₁ = 0` (∂Φ₁/∂x = 0). With the source's row convention, the linearized equation gives `(∂Φ₁/∂x)_new = 2 (∂Φ₁/∂x)_old`. Starting from a uniform Φ₁ hides this. It is the same pattern as numerical-audit §8.19 (NMC crystal model, group 02).
- **Evidence (T3, private instrumented 1C run):** the separator Φ₁ gradient and spread stay **exactly 0** for the whole run, because Φ₁ starts uniform and the row multiplies an exact zero. The defect is therefore latent: it has no effect on the archived results. It would matter for any start from a non-uniform Φ₁ (a restart, a charge after discharge, or a different initial condition). Fix in corrected mode: impose `(∂Φ₁/∂x)_new = 0`, or better, drop Φ₁ from the separator entirely.

## D-2. Separator porosity at the boundary-node faces

- **Where:** j = 1 row 1 (L1138–L1143); j = SEP_NODE rows 1 and 3, west side (L1239–L1256, L1337–L1354).
- **What:** these face fluxes lie in the separator but use `ε` (0.5) with separator transport coefficients. The separator nodes on the other side of the same faces (j = 2 and j = SEP_NODE−1) use `ε_sep` (0.39). The two nodes sharing a face therefore see fluxes that differ by a factor of ε_sep/ε = 0.78, so lithium and charge are not conserved across those faces.
- **Evidence (T3, private instrumented 1C run):**
  - Solid lithium gain equals I·t/F to 7 digits, and total lithium gain / (I·t/F) → 0.99974. Global lithium is conserved to within the salt loss below.
  - The mismatch acts at the two zero-volume nodes. The separator interior balances with ε_sep coefficients, while the boundary rows force the ε-based expression to carry the applied flux. In effect the separator *transports* with ε = 0.5 but *stores* with ε_sep = 0.39. Its steady concentration and potential gradients are therefore ε_sep/ε = 0.78 of the intended values: separator polarization is underestimated by about 22 %.
  - During start-up, while gradients form, the two nodes' creation and destruction don't cancel. The total anion inventory (which should be constant) drops 0.49 % in the first 18 s and settles at −0.520 % by about 170 s, then stays there to within 1 × 10⁻⁵.
- **Fix in corrected mode:** use ε_sep for every separator face. Check it with the anion-conservation test (T7), which should hold to machine precision.

## D-3. No voltage cutoff; the run ends on a NaN

- **Where:** main loop exits (L708–L724).
- **Evidence (L0, private oracle):** all five archived C-rate runs end with `EXIT BECAUSE delC ISNAN`, and the last output row has NaN voltage. Reruns reproduce this exactly.
- **Evidence (T3, private instrumented 1C run):** the first NaN appears at t = 3599 s in the reaction rate itself. Every cathode node has `c_s` slightly above `c_s,max` (θ = 1.00016 at the interface node), so `(c_s,max − c_s)^0.5` in the exchange current is NaN. With no solid diffusion and a flat OCP, lithiation is almost uniform, so the whole electrode reaches θ = 1 at once, at 100 % of theoretical capacity (3600 s at 1C). A single linearized step then overshoots θ = 1.
- **Fix in corrected mode:** voltage cutoff (the protocol engine, T9) and a guard on θ ∈ (0, 1) inside the kinetics.

## D-4. Precision and implicit typing

- Constants given as single-precision literals (see [parameters.md](parameters.md), "stored" column).
- `ex_1` in `echem_rxn` (L262) is implicitly `REAL`, so the exchange current density is rounded to single precision on every evaluation. This also affects the FD reaction derivatives.
- `last_write_time` in the main program is implicitly `INTEGER` (L683, L702), so the output spacing is 17 s rather than the intended 18 s. **Confirmed:** archived rows fall at 0, 18, 35, 52, 69 s …
- **Measured effect of `ex_1` (T3):** declaring `ex_1` double precision changes no digit of the archived 1C `Time_Voltage.txt` (5 printed decimals), so the effect is below 10⁻⁵ V. Faithful mode still reproduces the rounding so that raw values match the oracle closely.
- `10**(0.966)` in `rxn_k` (L98) is evaluated in single precision.
- `kappa` and `N_w_el` in `fillmat` are implicitly `INTEGER` but unused.

## D-5. Discarded Redlich–Kister sum in the OCP

- **Where:** L198–L201.
- **What:** computed on every OCP call (about 10 calls per node per step), then overwritten by the arctangent fit. It doesn't affect results but costs time, and it divides by `(2θ−1)^(1−k)`, which is singular at θ = 0.5.

## D-6. Finite-difference reaction derivatives

- **Where:** L1060–L1083.
- **What:** an absolute step of 10⁻⁶ in `c` (0.1 % of c⁰), in `c_s` (≈ 4 % of c_s,max, and 100 % of c_s at the start, where a one-sided difference is used), and in Φ (1 µV). Analytic derivatives are available.

## D-7. No Newton iteration within a time step

- **Where:** `bound_val_e` (L769–L827).
- **What:** one linear solve per step. The error from the nonlinearity is uncontrolled and depends on Δt.

## D-8. Inactive crystal-scale code

See [model.md §8](model.md#8-inactive-code-present-in-the-source-not-executed-in-any-archived-run). This has no effect on archived results because the code never runs.

## D-9. Volume fractions add up to more than 1

- `ε = 0.5` and `ε_AM = 0.8`. If `ε_AM` is meant as a fraction of the electrode volume, the solid and pore volumes sum to 1.3. It may instead be the active fraction *of the solid phase*. The equations use `ε_AM` directly in `a = 3ε_AM/R_p` and in the capacity. Needs the author's intent.

## D-10. Out-of-bounds read in the output routine

- **Where:** `write_all_voltage`, L406, L410, L414, L422.
- **Evidence (T3):** a `-fcheck=bounds` build stops at L414: `Index '-39' of dimension 2 of array 'cprev' below lower bound of 1`. The values read (`cs_NJ_2`, and `c0` at L422) are unused or overwritten before use (L464), so results are unaffected. It is still undefined behaviour and must not be ported. The same pattern appears in the NMC sources (numerical audit §8.36).

## D-11. Diffusion current missing from the ionic-current residual

- **Where:** row 3 of the separator interior (L1576–L1578), the separator/cathode interface (L1352–L1354) and the cathode interior (L1664–L1666).
- **What:** the Jacobian includes the diffusion part of the ionic current, `−εF(z₊D₊+z₋D₋)/τ ∂c/∂x` (coefficients `dW(3,1)`, `dE(3,1)`). The residual `smG(3)`, however, contains only the migration part (`dW(3,3)·dcdxW(3)`, `dE(3,3)·dcdxE(3)`). With one linear solve per step, the equation actually enforced is `i_mig(c_new) + [i_diff(c_new) − i_diff(c_old)] = source`. The diffusion current enters only through its change over one step, so in practice the model uses Ohm's law, `i₂ = −κ ∂Φ₂/∂x`, and drops the diffusion potential. The collector row (L1444) is the only one that includes the diffusion term.
- **Evidence (T3, private instrumented 1C run):** on the separator face between nodes 11 and 12, after start-up:
  - migration current = **0.780 I**, which is the ε_sep/ε factor of D-2, so this is the quantity the model conserves
  - diffusion current = −1.170 I
  - the full dilute-solution current is −0.39 I instead of I

  For a binary electrolyte with t₊ = 0.25 and zero anion flux, a consistent solution has migration = 2 I and diffusion = −I. The electrolyte potential gradient in the archived model is therefore about 2.6 times too small in the separator, and the concentration and potential fields are mutually inconsistent.
- **Fix in corrected mode:** include the full current in every residual row. Check it with a steady-state test (zero anion flux, so i₂ = F N₊ in the separator) and with `check_jacobian`, which cannot catch this defect alone because the Jacobian is right and the residual is wrong. A residual-consistency test is needed: G must equal −F(c) computed independently.
