# Changelog

## 0.4.0 (2026-10-03)

A crystal scale in corrected mode; faithful mode is unchanged, and the uniform model's output is byte-identical to 0.3.0 except in the two driver cases fixed below.

- **Crystal model** (`particle_model = 'crystal'`, docs/model.md §12): each cathode volume holds a crystal with solid-state diffusion, `crystal_shape` = `'sphere'` (default), `'cylinder'` or `'slab'`, solved together with the electrode by a condensed Newton step. The original's inactive crystal code is not ported (D-8); the α/β phase change is not included.
- **Inputs:** `particle_model`, `crystal_shape` and `D_c` (default 8 × 10⁻¹⁴ cm²/s, the original's `diff_c`) in `&active`; `nj_crystal` (default 21) in `&cell`.
- **Notebook** `04_crystal_scale`: profiles inside the crystals, capacity against rate and D_c, the uniform limit and the three shapes.
- Requires bandsolver 0.1.2.
- **Driver fixes** (corrected mode, all three languages): a cc discharge ends only at its `Vmin` and a charge only at its `Vmax`, as documented (both bounds were applied, so a step starting beyond the other bound stopped at once as a cutoff); when a time step cannot be solved, the exit row is the last converged sub-step, not the state at the start of the time step (the capacity of the partial step was lost). Found while porting the Zn/MnO₂ model.

## 0.3.0 (2026-10-03)

Corrected mode is reformulated so that physical limits are reached smoothly; faithful mode is unchanged (still byte-identical to the original).

- **Corrected mode, log variables (model.md §10):** the unknowns are u = ln(c/c_bulk), Φ₁, Φ₂ and the particles' log-odds s = ln(θ/(1−θ)), so c > 0 and 0 < θ < 1 by construction. Ion fluxes use exponential fitting (Scharfetter–Gummel). The electrolyte can run out and particles can fill or empty without clipping. Salt and lithium are conserved to round-off. This supersedes the exchange-current regularization D-13.
- **Corrected mode, OCP tails (D-16):** the LFP open-circuit potential has ideal-solution tails at empty and full (θ_e = 10⁻⁴), so particles approach either end asymptotically. They change U by at most 0.25 mV between 1 % and 99 % lithiation. Discharges change by at most 0.02 mV away from the ends and 0.55 mV on the final drop to the cutoff (1C). After a full charge the rest voltage relaxes to 3.76 V instead of 3.43 V.
- **Input:** new `kappa_bg` in `&electrolyte` (default 10⁻⁸ S/cm), the solvent's background conductivity, which keeps Φ₂ defined where the salt is exhausted.
- **Driver:** sub-steps halve down to 10⁻¹⁰ s and double again after each success; a time step gives up after 200 failures. When a step cannot be solved, the exit reason names the limit reached (`electrolyte_depleted`, `particles_full`, `particles_empty`), or `solver_fail` if none applies, and the last row reports the state at the start of that step.
- Python: `lfp_model.logcore` and `lfp_model.logmodel.LogModel` (with `conc` and `cs` to convert a state to concentrations); `Result.final_state` in corrected mode now holds (u, Φ₁, Φ₂, s). The package's `__version__` was stale at 0.1.0 and is now 0.3.0.

## 0.2.0 (2026-09-28)

Changes to inputs and to corrected-mode results; faithful mode is unchanged (still byte-identical to the original).

- **Input:** the cathode thickness is now `L_cath_um` in µm (default 24) in the `&cell` group, replacing `L_cath` in cm. The original computes the thickness as `24 * 1.0d-4`, which a decimal literal in cm cannot reproduce exactly. Replace `L_cath = 24.0d-4` with `L_cath_um = 24` in existing input files.
- **Corrected mode, voltage reference (D-15):** potentials are referenced to the lithium foil (0 V). The cell voltage now includes the foil's Nernst term, V = Φ₁(collector) − U_Li − η_Li. Constant-voltage steps hold this voltage. The output has a new tenth column, `Li_Nernst` (mV).
- **Corrected mode, OCP (D-14):** the LFP open-circuit potential includes the electrolyte Nernst term (RT/F)·ln(c/c_bulk), with its exact derivative in the Jacobian.
- Together these shift a corrected 1C discharge by −1.6 mV and a 2C discharge by −3.2 mV. Capacities are unchanged.
- Python: `kinetics.ocp(p, cs, c=None)`, `kinetics.li_foil`, `simulate.foil_shift` and `simulate.foil_referenced`.

## 0.1.0 (2026-09-28)

First public version.

- Fortran (`fortran/lfp.f90`), C++ (`cpp/lfp.cpp`) and Python (`lfp_model`) implementations that read one namelist input file and write the same output.
- Faithful mode reproduces the original research program's output byte for byte.
- Corrected mode fixes the defects listed in `docs/deviations.md`:
  - diffusion current in the ionic-current residual
  - separator porosity on the boundary faces
  - voltage cutoffs
  - Newton iteration with exact kinetic derivatives
  - the Li-foil boundary-row sign
  - double precision throughout
  - active fraction of the solid phase
  - Butler–Volmer lithium counter electrode
- Cycling protocols in corrected mode: constant current, constant voltage and rest steps, repeated `cycles` times (`docs/protocol.md`).
- Test suite (conservation, independent residual, Jacobian, convergence, protocols, cross-language agreement), CI, and three tutorial notebooks.
