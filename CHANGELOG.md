# Changelog

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
