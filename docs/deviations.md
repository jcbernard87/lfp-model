# Deviations from the original research code

The faithful mode (`mode = "faithful"`) reproduces the original program, including every item below. The corrected mode (`mode = "corrected"`, the default once T8 is done) applies the fixes marked **fixed**. Each fix is a separate commit that states the effect on results.

Status values: **candidate** (suspected from reading the source), **confirmed** (demonstrated by a test or run), **fixed**, **kept** (reviewed and left as is, with the reason).

| ID | Status | Summary |
|---|---|---|
| D-1 | candidate | Li-foil face (j = 1): the solid-potential row doubles the old gradient instead of zeroing it |
| D-2 | candidate | Separator faces next to the boundary nodes use the cathode porosity `ε`, so flux doesn't balance across faces 1|2 and (SEP−1)|SEP |
| D-3 | confirmed | No voltage cutoff: every discharge ends on a NaN |
| D-4 | candidate | Precision and typing: single-precision literals, implicitly single-precision `ex_1`, implicitly integer `last_write_time` |
| D-5 | candidate | The OCP routine computes and discards a 31-term Redlich–Kister sum on every call |
| D-6 | candidate | Reaction derivatives use absolute FD steps of 10⁻⁶ (≈ 4 % of c_s,max, 100 % of the initial c_s) |
| D-7 | candidate | One linearized solve per time step (no Newton iteration); nonlinearity error is not controlled |
| D-8 | candidate | The inactive crystal-scale code has undeclared and misspelled variables and inconsistent constants |
| D-9 | candidate | Porosity 0.5 plus active-material fraction 0.8 add up to more than 1 |

## D-1. Li-foil face solid-potential row sign

- **Where:** `fillmat`, j = 1, row 2 (L1168–L1169).
- **What:** the row is intended to impose `i₁ = 0` (∂Φ₁/∂x = 0). With the source's row convention, the linearized equation gives `(∂Φ₁/∂x)_new = 2 (∂Φ₁/∂x)_old`. Starting from a uniform Φ₁ hides this. It is the same pattern as numerical-audit §8.19 (NMC crystal model, group 02).
- **Evidence needed (T3):** separator Φ₁ gradient over a run; whether roundoff grows.

## D-2. Separator porosity at the boundary-node faces

- **Where:** j = 1 row 1 (L1138–L1143); j = SEP_NODE rows 1 and 3, west side (L1239–L1256, L1337–L1354).
- **What:** these face fluxes lie in the separator but use `ε` (0.5) with separator transport coefficients. The separator nodes on the other side of the same faces (j = 2 and j = SEP_NODE−1) use `ε_sep` (0.39). The two nodes sharing a face therefore see fluxes that differ by a factor of ε_sep/ε = 0.78, so lithium and charge are not conserved across those faces.
- **Evidence needed (T3):** a per-step lithium balance (electrolyte + solid) against the charge passed.

## D-3. No voltage cutoff; the run ends on a NaN

- **Where:** main loop exits (L708–L724).
- **Evidence (L0, private oracle):** all five archived C-rate runs end with `EXIT BECAUSE delC ISNAN`, and the last output row has NaN voltage. Reruns reproduce this exactly.
- **Evidence needed (T3):** which expression produces the first NaN (suspected: `(c_s,max − c_s)^0.5` in the exchange current once `c_s` or `c_s + step` exceeds `c_s,max`).

## D-4. Precision and implicit typing

- Constants given as single-precision literals (see [parameters.md](parameters.md), "stored" column).
- `ex_1` in `echem_rxn` (L262) is implicitly `REAL`, so the exchange current density is rounded to single precision on every evaluation. This also affects the FD reaction derivatives.
- `last_write_time` in the main program is implicitly `INTEGER` (L683, L702), so the output spacing is 17 s rather than the intended 18 s.
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
