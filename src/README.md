# `src/` — Standalone Fortran modules

This directory contains the standalone EOS modules. They compile and
run independently of PFLOTRAN, using only standard Fortran 2008 and a
single Makefile.

## Build

```bash
make           # build all six test executables
make test      # run PR + Peneloux EOS anchor check (NIST data)
make sweep     # run 440-point P-T sweep, write air_eos_sweep.csv
make crosscheck   # Python independent verification of the sweep
make henry        # Henry's law test
make watervapour  # water-vapour partitioning test
make properties   # viscosity + ideal-gas enthalpy test
make wrapper      # PFLOTRAN-facing wrapper test
make clean        # remove build artifacts
```

All test executables should print uniform "PASS" output and complete
in well under a second each.

## Module map

| Module | Purpose | Key references |
|---|---|---|
| `air_eos_pr_module.F90` | Peng-Robinson + Peneloux density, fugacity, h_dep | Peng & Robinson 1976; Péneloux et al. 1982; Lemmon et al. 2000 |
| `air_henry_module.F90` | Henry's law O₂/N₂/CO₂/Ar with salting-out | Sander 2015; Weisenberger & Schumpe 1996 |
| `air_watervapour_module.F90` | IAPWS-IF97 + Poynting + brine activity | Wagner et al. 2000 |
| `air_properties_module.F90` | Lemmon-Jacobsen viscosity + Shomate h_ig | Lemmon & Jacobsen 2004; NIST WebBook |
| `air_eos_module.F90` | Unified PFLOTRAN-facing wrapper | This work |

Each module has a matching `test_<name>.F90` driver in this directory.

## Cross-validation

`cross_validate.py` is an independent Python implementation of the PR
+ Peneloux EOS, used to verify the Fortran output to bit precision.

```bash
make sweep         # writes air_eos_sweep.csv from Fortran
python3 cross_validate.py air_eos_sweep.csv
```

## Integration

To use these modules in PFLOTRAN, see [`../pflotran_integration/`](../pflotran_integration/)
and [`../docs/WORKFLOW.md`](../docs/WORKFLOW.md).
