# LFP model parameters

These are the default values in the original research code. The source line (`L…`) refers to the private reference copy (see [model.md](model.md)).

**About the "stored" column.** Several constants are written as default-kind (single-precision) Fortran literals, such as `8.314` rather than `8.314d0`, and are then assigned to double-precision variables. The value the original program actually uses is the single-precision rounding, shown under "stored" where it differs. The faithful mode reproduces these stored values. The corrected mode uses the intended decimal values. **[D-4]**

## Cell geometry and mesh

| Name | Symbol | Value | Unit | Source | Notes |
|---|---|---|---|---|---|
| `xmax` | L_cath | 24 × 10⁻⁴ | cm | L124 | cathode thickness, 24 µm |
| `len_sep` | L_sep | 25 × 10⁻⁴ | cm | L109 | separator thickness, 25 µm (the source comment asks whether µm or cm) |
| `NJ` | | 101 | | L35 | total nodes |
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
| `diff_c` | D_s | 8.0 × 10⁻¹⁴ | cm²/s | L97 | used only by the inactive crystal scale |

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
| anode kinetic constant | 10⁻⁶ | | L241 | Li-foil exchange current, output only |
