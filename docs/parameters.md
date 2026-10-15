# LFP model parameters

These are the default values in the original research code. The source line (`L…`) refers to the private reference copy (see [model.md](model.md)).

**About the "stored" column.** Several constants are written as default-kind (single-precision) Fortran literals, such as `8.314` rather than `8.314d0`, and are then assigned to double-precision variables. The value the original program actually uses is the single-precision rounding, shown under "stored" where it differs. The faithful mode reproduces these stored values. The corrected mode uses the intended decimal values. **[D-4]**

## Cell geometry and mesh

| Name | Symbol | Value | Unit | Source | Notes |
|---|---|---|---|---|---|
| `xmax` | L_cath (input `L_cath_um`, µm) | 24 × 10⁻⁴ | cm | L124 | cathode thickness, 24 µm. The original computes `24 * 1.0d-4`, which differs from the literal 24.0 × 10⁻⁴ in the last bit, so the input is the thickness in µm and L_cath = L_cath_um × 10⁻⁴ cm |
| `len_sep` | L_sep | 25 × 10⁻⁴ | cm | L109 | separator thickness, 25 µm (the source comment asks whether µm or cm) |
| `NJ` | | 101 | | L35 | total nodes |
| (new) `nj_crystal` | | 21 | | — | crystal model: nodes across a crystal, center and surface included (the original's `NJ_c` = 51 was never used). 21 nodes give the 2C capacity within 3 × 10⁻⁴ equivalents of 81 nodes (validation.md §6) |
| (new) `nja` | | 21 | | — | agglomerate model (`&agglomerate`): nodes along an agglomerate radius, center and surface included (zero-volume) |
| (new) `R_agg` | | 5 × 10⁻⁴ | cm | — | agglomerate model: agglomerate radius (5 µm, a typical spray-dried LiFePO₄ secondary particle); the crystals inside have radius `R_p` |
| (new) `eps_agg` | | 0.2 | | — | agglomerate model: agglomerate porosity; ε + f_AM(1−ε)/(1−ε_agg) ≤ 1 |
| (new) `tau_agg`, `sigma_agg` | | −1, −1 | , S/cm | — | agglomerate model: pore tortuosity and solid conductivity inside the agglomerates; ≤ 0 means ε_agg^−½ and `sigma` |
| `SEP_NODE` | | 22 | | L108 | separator/cathode interface node |
| `eps` | ε | 0.5 | | L126 | cathode porosity |
| `volfrac_AM` | ε_AM (`eps_AM`) | 0.8 (stored 0.800000011920929) | | L125 | faithful mode: active-material volume fraction, as in the original (ε + ε_AM = 1.3 > 1, **[D-9]**) |
| (new) | f_AM (`f_AM`) | 0.8 | | — | corrected mode: active fraction of the solid phase; active volume fraction f_AM·(1−ε) = 0.4 **[D-9]** |
| `tortuosity` | τ | ε^(−1/2) | | L128 | Bruggeman |
| `eps_sep` | ε_sep | 0.39 (stored 0.38999998569488525) | | L107 | |
| `sep_tortuosity` | τ_sep | 4.0 | | L110 | |

## Electrolyte

| Name | Symbol | Value | Unit | Source | Notes |
|---|---|---|---|---|---|
| `diff` | D | 2.00 × 10⁻⁶ | cm²/s | L57 | the only difference from `LFP_03042021.f95`, which uses 3.5 × 10⁻⁶ |
| `transference_num_cat` | t₊ | 0.25 | | L58 | |
| `cbulk` | c⁰ | 0.001 (stored 0.0010000000474974513) | mol/cm³ | L61 | 1 M |
| `z_cat`, `z_an` | z₊, z₋ | +1, −1 | | L63–L64 | |
| (new) `kappa_bg` | κ_bg | 10⁻⁸ | S/cm | — | corrected mode (v0.3.0): the solvent's own ionic conductivity (neat carbonates are about 10⁻⁹ to 10⁻⁷ S/cm), about 10⁻⁶ of the 1 M electrolyte's; keeps Φ₂ defined where the salt is exhausted (model.md §10) |

## Active material (LiFePO₄)

| Name | Symbol | Value | Unit | Source | Notes |
|---|---|---|---|---|---|
| `sigma` | σ | 3.0 × 10⁻³ | S/cm | L70 | |
| `molar_mass_AM` | M | 125.759 (stored 125.75900268554688) | g/mol | L71 | |
| `density_AM` | ρ | 3.6 (stored 3.5999999046325684) | g/cm³ | L72 | |
| `Qmth` | Q_th | 0.170 (stored 0.17000000178813934) | Ah/g | L73 | |
| derived `cimax` | c_s,max | ρ Q_th · 3600/F ≈ 0.02283 | mol/cm³ | L233 | |
| `xmax_c` | R_p | 200 × 10⁻⁷ | cm | L96 | particle radius (200 nm); sets `a` |
| derived `spec_a` | a | 3 ε_AM/R_p = 1.2 × 10⁵ | cm⁻¹ | L127 | |
| `rxn_k` | k | 10⁻⁸ · 10^0.966 (stored 10⁻⁸ × 9.24698257446289) | mol^(−1/2)·cm^(5/2)/s | L98 | 10^0.966 is evaluated in single precision, with the exponent itself rounded to float32 first; confirmed by the byte-exact faithful reproduction |
| `alpha_a`, `alpha_c` | α_a, α_c | 0.5, 0.5 | | L148 | |
| `diff_c` | D_c (`D_c`) | 8.0 × 10⁻¹⁴ | cm²/s | L97 | used only by the inactive crystal scale of the original (its sweep substitution is commented out); in v0.4.0 the default solid diffusivity of the crystal model (model.md §12) |
| (new) | `particle_model` | `'uniform'` | | — | corrected mode: `'uniform'` particles (§10) or `'crystal'` (solid diffusion, model.md §12) |
| (new) | `crystal_shape` | `'sphere'` | | — | crystal model: `'sphere'`, `'cylinder'` or `'slab'`; a = (k+1)ε_AM/R_p with k = 2, 1, 0 |

## Constants

| Name | Value | Stored | Source |
|---|---|---|---|
| `Rigc` | 8.314 J/(mol·K) | 8.314000129699707 | L50 |
| `Temp` | 298 K | exact | L51 |
| `Fconst` | 96485 C/mol | exact | L52 |
| `PI` | 3.141592654 | single-precision `REAL` (unused) | L53 |

## Operating conditions and numerics

| Name | Value | Unit | Source | Notes |
|---|---|---|---|---|
| `C_rate` | 0.1, 0.2, 0.5, 1, 2 in the archived runs | 1/h | L116 | set by the run generator as a single-precision literal, so 0.1 and 0.2 are stored as 0.10000000149… and 0.20000000298… |
| derived `c_specific` | Q_th × C_rate | A/g | L130 | |
| derived `c_density` | I = c_specific · L_cath · ε_AM · ρ | A/cm² | L131 | 1.175 × 10⁻³ A/cm² at 1C |
| `Phi_1_init` | 3.6 (stored 3.5999999046325684) | V | L79 | |
| `Phi_2_init` | 0 | V | L80 | |
| `c0_init` | 0.001 | mol/cm³ | L81 | unused; the initial electrolyte concentration is `cbulk` |
| `cs_init` | 1.0 × 10⁻⁵ | mol/cm³ | L82 | θ₀ ≈ 4.4 × 10⁻⁴ (fully delithiated) |
| `tmax`, `Numbertimesteps` | 36 000 s, 36 000 | | L35 | Δt = 1 s |
| FD reaction-derivative steps | 10⁻⁶ (absolute) | | L1061–L1064 | **[D-6]** |
| (new) `dt_mode` | `'fixed'` | | — | corrected mode (#18): `'adaptive'` chooses the time step from the local error estimate (docs/protocol.md) |
| (new) `rtol`, `atol` | 10⁻⁴, 10⁻⁴ | | — | adaptive: local error tolerances on the state (ln c, potentials in V, log-odds of θ) |
| (new) `dt_min`, `dt_max` | 10⁻⁶, 0 | s | — | adaptive: step bounds; `dt_max = 0` means `write_interval` |
| (new) `restart_read`, `restart_write` | `''` | | — | `&output`, corrected mode: start from, or write the final state to, a restart file (docs/protocol.md; #23) |
| (new) `columns` | `''` | | — | `&output`, corrected mode: the time-series columns to write; `''` = all (docs/protocol.md; #22) |
| (new) `profile_file`, `profile_vars` | `''` | | — | `&output`, corrected mode: write profiles against x to this file; their variables (`''` = all) |
| (new) `profile_every`, `profile_times` | 0, `''` | , s | — | `&output`: a profile with every N-th output row; at the first row at or after each listed time |
| anode kinetic constant | 10⁻⁶ | | L241 | Li-foil exchange current, output only |
