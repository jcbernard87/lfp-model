# lfp-model

A one-dimensional model of a LiFePO₄ half cell (lithium foil | separator | porous LiFePO₄ cathode | current collector). It is solved with [bandsolver](https://github.com/jcbernard87/bandsolver), an implementation of Newman's BAND method.

The model uses dilute-solution transport of a binary electrolyte, electronic conduction in the cathode, Butler–Volmer insertion kinetics, and uniform-concentration particles. It is discretized with finite volumes and advanced by backward Euler, with Newton's method at every step. It runs constant-current, constant-voltage and rest protocols, including cycling.

It comes in three implementations that read the same input file and write the same output:

| Implementation | Where | Solver | Use it for |
|---|---|---|---|
| Fortran program | [`fortran/lfp.f90`](fortran/lfp.f90) (one file) | links bandsolver | batch runs |
| C++ program | [`cpp/lfp.cpp`](cpp/lfp.cpp) (one file) | its own embedded copy of bandsolver's C++ core; no dependencies | batch runs, embedding |
| Python package | [`python/lfp_model`](python/lfp_model), [`notebooks/`](notebooks) | bandsolver (Python) | examples, quick studies, notebooks |

This repository shares the **model only**. It contains no simulation results or experimental data; run the model to generate results.

## Quick start

**C++** (only a C++17 compiler is needed):

```sh
c++ -std=c++17 -O2 -ffp-contract=off cpp/lfp.cpp -o lfp_cpp
./lfp_cpp input/default.nml          # writes Time_Voltage.txt
```

**Fortran and C++ with CMake** (bandsolver is downloaded automatically):

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
build/fortran/lfp_f input/default.nml
```

With a local bandsolver checkout, add `-DFETCHCONTENT_SOURCE_DIR_BANDSOLVER=/path/to/bandsolver`.

**Python:**

```sh
pip install "bandsolver @ git+https://github.com/jcbernard87/bandsolver@v0.1.1"
pip install -e ".[notebooks]"
python -m lfp_model input/default.nml
```

```python
from lfp_model.params import Params
from lfp_model.simulate import run

r = run(Params(mode="corrected", C_rate=1.0))
t_h, volts = r.array[:, 0], r.array[:, 1]
```

## Input

One Fortran-namelist file, [`input/default.nml`](input/default.nml), holds every parameter: geometry, electrolyte, active material, operation, numerics and output. [docs/parameters.md](docs/parameters.md) lists each parameter with its unit and origin. The Fortran program reads the file natively; the C++ and Python implementations parse the same format.

A cycling protocol goes in a `&protocol` group ([docs/protocol.md](docs/protocol.md)):

```fortran
&protocol
  steps  = 'cc C=1 Vmin=2.5; rest t=1800; cc C=-1 Vmax=4.0; cv V=4.0 Imin=0.05; rest t=1800'
  cycles = 2
/
```

## Two modes

The model is a port of the author's PhD research code. The port was checked for errors while it was being made:

- **`mode = 'corrected'`** (recommended) fixes every defect found. The fixes are listed with evidence and their measured effect in [docs/deviations.md](docs/deviations.md).
- **`mode = 'faithful'`** reproduces the original program exactly, defects included. It is kept for comparison with earlier work. In faithful mode all three implementations reproduce the original program's output files byte for byte ([docs/validation.md](docs/validation.md)). Faithful mode runs only the original constant-current discharge.

`default.nml` selects faithful mode so that a first run reproduces the original. Set `mode = 'corrected'` for new work.

## Documentation

- [docs/model.md](docs/model.md): equations, boundary conditions, discretization and time stepping
- [docs/parameters.md](docs/parameters.md): every parameter, with value, unit and origin
- [docs/protocol.md](docs/protocol.md): protocol steps and the output columns
- [docs/deviations.md](docs/deviations.md): defects found in the original code, with evidence and fixes
- [docs/validation.md](docs/validation.md): reproduction, cross-language agreement, conservation and convergence tests, with measured numbers
- [notebooks/](notebooks): quickstart, the model explained step by step, and rate capability with a full cycle

## Tests

```sh
cmake --build build && python -m pytest
```

The public tests compute everything they need. They cover cross-language agreement, an independent residual, the Jacobian, conservation, equilibrium, time and mesh convergence, and protocols. The tests in `private_checks/` compare against the original program and its archived runs. They run only when `LFP_ORACLE_DIR` points at that material, which is not distributed.

## Status

Version 0.1.0, in development. The inactive crystal-scale (solid-diffusion) submodel of the original code is not ported (deviation D-8). Deviations D-9 and D-12 involve modelling choices awaiting the author's confirmation.

## License and citation

BSD 3-Clause; see [LICENSE](LICENSE). If you use this model, please cite [bandsolver](https://github.com/jcbernard87/bandsolver) and Newman's method: J. Newman, *Ind. Eng. Chem. Fundam.* 7, 514 (1968); J. Newman and K. E. Thomas-Alyea, *Electrochemical Systems*, 3rd ed., Appendix C.
