# lfp-model

A LiFePO₄ porous-electrode battery model with crystal-scale lithium transport and a phase-change model, solved with [bandsolver](https://github.com/jcbernard87/bandsolver) (Newman's BAND method).

**Status: in development.** The model is being ported from the author's PhD research code. Work in progress is tracked in `LFP_LOOP.md`.

The repository contains three implementations that read the same input file and write the same outputs:

| Implementation | Location | Solver |
|---|---|---|
| Fortran program | `fortran/lfp.f90` | links bandsolver |
| C++ program | `cpp/lfp.cpp` | embedded copy of bandsolver's C++ core |
| Python package + notebooks | `python/lfp_model`, `notebooks/` | bandsolver |

This repository shares the model only. It does not include simulation results or experimental data; run the model to generate results.

## License

BSD 3-Clause; see [LICENSE](LICENSE).
