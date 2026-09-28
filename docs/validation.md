# Validation

All numbers below were measured before deviation D-9 was applied to corrected mode, except the test results, which are re-run on every change. D-9 halves the active volume fraction, so per-area quantities (currents, a) change, but none of the checks depend on it. All numbers were measured on macOS arm64 (gfortran 14.2, Apple clang 16, Python 3.13, numpy 2.5, bandsolver 0.1.1) on 2026-09-27. The public checks run in CI (`pytest`) and recompute everything; no result files are stored in this repository.

## 1. Reproduction of the original research code (private)

The original program and its archived outputs are kept privately (they contain third-party-derived solver code and research results). With `LFP_ORACLE_DIR` pointing at that material, `private_checks/` compares the faithful mode of every implementation with the archived `Time_Voltage.txt` files of five constant-current discharges (0.1C, 0.2C, 0.5C, 1C, 2C):

| Implementation | Result | Time at 1C |
|---|---|---|
| original program, rebuilt (`gfortran -O0`) | byte-identical to the archived files | 1.0 s |
| Python (`lfp_model`, faithful) | byte-identical | 3.4 s |
| Fortran (`fortran/lfp.f90`, faithful) | byte-identical | 0.38 s |
| C++ (`cpp/lfp.cpp`, faithful) | byte-identical | 0.23 s |

Byte-identical output required reproducing the original's single-precision literals, its implicit single-precision exchange current, its integer output timer, and building without fused multiply-add contraction (`-ffp-contract=off`). See [deviations.md](deviations.md) D-4.

## 2. Cross-language agreement (public, `tests/test_cross_language.py`)

| Comparison | Mode | Result |
|---|---|---|
| Fortran vs Python, 2C and 0.5C | faithful | identical files |
| C++ vs Python, 2C and 0.5C | faithful | identical files |
| Fortran vs Python, 2C and 0.5C | corrected | agree to rtol 1e-5 (rounding of the printed values) |
| C++ vs Python, 2C and 0.5C | corrected | agree to rtol 1e-5 |
| C++ vs Fortran, 2C and 0.5C | corrected | identical files |
| Full protocol cycle: Fortran, C++, Python | corrected | Fortran and C++ files identical; Python agrees to rtol 1e-5 |

## 3. Physics and numerics (public, `tests/test_physics.py`)

| Check | Corrected mode | Faithful mode |
|---|---|---|
| Assembled residual vs an independent face-by-face implementation of [model.md](model.md) (perturbed mid-discharge state) | agree to 3 × 10⁻¹⁶ relative | cation and current rows differ by O(1): D-2, D-11 |
| `bandsolver.check_jacobian` at a mid-discharge state | max error 3.9 × 10⁻⁶ (finite-difference noise) | sign error at node 1, solid-potential row: D-1 |
| Salt inventory over 1800 s at 1C | constant to 4 × 10⁻¹⁶ | −0.520 %: D-2 |
| Lithium into the solid / (I t/F) | 1 − 1 × 10⁻¹⁶ | 1 − 3 × 10⁻¹⁵ |
| Ionic current across the separator | = I to 10⁻⁹ | −0.39 I (D-11) |
| Rest at equilibrium (no current, Φ₁ = U) | no change to 10⁻¹² | not tested |
| Time-step convergence (Li-face concentration at t = 8 s, 2C, Δt = 1 … 0.0625 s) | orders 0.96, 0.98, 0.99 (backward Euler: 1) | |
| Mesh convergence (cell voltage at 900 s, 2C, 26 … 401 nodes) | orders 2.11, 2.05, 2.02 | |
| End of discharge | `cutoff_low`, located within 0.1 mV of 2.5 V; 0.7964–0.7973 electron equivalents (2C, 1C) | ends on a NaN (D-3) |
| Protocol cycle: 2C discharge, rest, 1C charge to 4.0 V, CV to 0.05C, rest (`tests/test_protocol.py`) | completes; the charge returns exactly the initial intercalated lithium (−θ₀·x_max, within 2 %); the final rest relaxes to the OCP within 2 mV | not applicable |

## 4. Effect of each fix on results

Measured by reverting one fix at a time in corrected mode, for a 1C discharge. ΔV is the change in cell voltage over the first 0.9 h, relative to corrected mode. These were measured before D-12 was introduced. D-12 changes only the reported counter-electrode overpotential, which is identical in every row, so the relative effects stand. D-12 itself shifts the cell voltage by −24.5 mV (0.1C) to −82 mV (2C); see [deviations.md](deviations.md#d-12-lithium-counter-electrode-overpotential).

| Reverted | Mean ΔV | Max \|ΔV\| | Notes |
|---|---|---|---|
| D-1 (Li-face Φ₁ row sign) | n/a | n/a | With Newton iteration the doubled-gradient row makes the solve fail at once. The defect was latent in the original only because it does one solve per step from an exactly uniform Φ₁. |
| D-2 (ε on the separator boundary faces) | +0.37 mV | 0.40 mV | |
| D-11 (diffusion current missing from the residual) | +1.90 mV | 1.94 mV | the largest single effect |
| D-6 (finite-difference reaction derivatives) | 0.00 mV | 0.00 mV | Changes only Newton convergence. Near full lithiation, Newton then fails before the cutoff (0.7937 vs 0.7973 equivalents). |
| all (faithful mode: every defect, one solve per step) | +1.85 mV | 2.37 mV | |

For this cell (24 µm cathode, 1 M electrolyte) electrolyte transport is not limiting, so the transport defects shift the voltage by only millivolts. They would matter more for thick electrodes, high rates or dilute electrolytes.

## 5. Not covered yet

- Comparison with experimental data. This is done privately and never published with the repository.
- The crystal-scale (solid-diffusion) model, which the original never ran (D-8).
