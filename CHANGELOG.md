# Changelog

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
