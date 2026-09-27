# Cycling protocols

In corrected mode the cell is driven by a **protocol**: a list of steps that run in order, optionally repeated. Faithful mode ignores the protocol and runs the original constant-current discharge.

The protocol is given in the `&protocol` group of the input file, read identically by the Fortran, C++ and Python implementations:

```fortran
&protocol
  steps  = 'cc C=0.5 Vmin=2.5; rest t=1800; cc C=-0.5 Vmax=4.0; cv V=4.0 Imin=0.05; rest t=1800'
  cycles = 3
/
```

Without a `&protocol` group (or with `steps = ''`) the protocol is a single step, `cc C=<C_rate> Vmin=<V_min> Vmax=<V_max>`, i.e. the classic discharge.

## Steps

Steps are separated by `;`. Each step is a keyword followed by `key=value` settings (case-insensitive, any order):

| Step | Settings | Ends when |
|---|---|---|
| `cc` | `C` = C-rate (1/h), **positive = discharge, negative = charge** (required); `t` = maximum duration [s]; `Vmin`, `Vmax` [V] (default: the global `V_min`, `V_max`) | the voltage reaches `Vmin` (discharge) or `Vmax` (charge), or `t` has elapsed |
| `cv` | `V` = held voltage [V] (required); `t` = maximum duration [s]; `Imin` = current magnitude, as a C-rate, below which the step ends | \|I\| ≤ `Imin`, or `t` has elapsed |
| `rest` | `t` = duration [s] (required) | `t` has elapsed |

A step with no possible end (for example `cv` with neither `t` nor `Imin`) is rejected. `cycles = n` runs the whole list n times.

**C-rate.** 1C is the current that would pass the theoretical capacity of the cathode in one hour: I = Q_th · L_cath · ε_AM · ρ (A/cm²).

## How steps are solved

- **Time step.** Every step uses Δt = `t_max/n_steps` (1 s by default). A step's final time step is shortened so the step ends exactly at `t`, and a voltage end condition is located inside a time step by sub-stepping (the same mechanism as the corrected-mode cutoff; the crossing is resolved to within one sub-step).
- **Constant current and rest** use the applied current directly (rest: I = 0).
- **Constant voltage.** The current is unknown: each time step finds the current I for which the cell voltage equals `V`, by secant iteration on I around the full Newton solve (starting from the previous step's current; |V − V_set| ≤ 10⁻⁹ V).
- **Cell voltage.** V = Φ₁(collector) − η_Li(I). The lithium counter electrode is treated as a symmetric Butler–Volmer interface, η_Li = (RT/(αF))·asinh(I/(2 i₀,Li)) with α = 0.5 and i₀,Li from the electrolyte concentration at x = 0. See deviation D-12: the original used (RT/F)·ln(I/i₀,Li), which is singular at zero current.

## Output

In corrected mode `Time_Voltage.txt` has the original seven columns plus:

| Column | Header | Meaning |
|---|---|---|
| 8 | Current (mA/cm2) | applied current density, positive on discharge |
| 9 | Step | 1-based index of the step in the expanded list (cycles × steps) |

Rows are written every `write_interval` seconds (default 18 s) and at the start and end of every step. The State column is `D` (I > 0), `C` (I < 0) or `R` (I = 0).

The run's exit reason is `end_of_protocol` when every step completed, or `nan`, `solver_fail` or `max_time` (99 h) otherwise. A single-step protocol reports how that step ended (`cutoff_low`, `cutoff_high`, `duration` or `current_limit`), which keeps the classic discharge's `cutoff_low`.
