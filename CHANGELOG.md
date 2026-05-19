# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-05-19

Initial release.

### Added

- Five standalone Fortran modules implementing the high-pressure air EOS:
  - `air_eos_pr_module.F90`: Peng-Robinson + Peneloux density,
    fugacity coefficient, enthalpy departure, with analytic derivatives.
    Calibrated against NIST air data (440-point P-T grid, max error ~10⁻⁷).
  - `air_henry_module.F90`: Henry's law for O₂/N₂/CO₂/Ar with Sander
    (2015) reference constants, van't Hoff temperature dependence,
    Setschenow/Weisenberger-Schumpe salting-out, and fugacity correction.
  - `air_watervapour_module.F90`: IAPWS-IF97 Region 4 saturation
    pressure, ideal-mixing water activity for NaCl brines, Poynting
    correction, coupled to bulk-gas fugacity coefficient.
  - `air_properties_module.F90`: Lemmon-Jacobsen (2004) viscosity
    (Chapman-Enskog dilute + 5-term residual polynomial) and NIST
    Shomate ideal-gas heat capacity / reference enthalpy.
  - `air_eos_module.F90`: PFLOTRAN-facing wrapper consolidating the
    above into a unified interface.

- Test driver suite (six executables) with NIST/literature anchors,
  finite-difference derivative checks, and a 440-point P-T sweep.

- Python cross-validation script reproducing the Fortran PR+Peneloux
  output to bit precision.

- PFLOTRAN integration:
  - Adapter routines in `eos_gas.F90` matching the function-pointer
    Dummy interfaces (~150 lines).
  - Setter routine `EOSGasSetPRPenelouxAir()` binding the pointers.
  - New input deck keyword `PR_PENELOUX_AIR` in `eos.F90` dispatcher.

- Regression test decks at increasing physical complexity:
  - Level 1: single-cell static EOS dispatch check.
  - Level 2: single-cell pressure transient.
  - Level 3: 1D column gas injection.
  - Level 5: 2D axisymmetric radial geometry (Yang et al. 2024 base case).

- Documentation:
  - `docs/WORKFLOW.md`: complete reproducible workflow.
  - `docs/INTEGRATION.md`: PFLOTRAN-specific integration steps.
  - `docs/methods.md`: publication-ready methods description.

### Verified

- Standalone EOS values agree with NIST sweep to better than 10⁻⁶
  across 0–200 bar, 0–100 °C.
- PFLOTRAN-integrated density matches the standalone reference to 0.03%
  at 80 bar, 40 °C.
- PFLOTRAN-integrated viscosity matches to numerical roundoff.
- Newton solver converges in ~2 iterations per timestep under flow.
- 1D injection simulation produces physically reasonable gas-front
  velocity (~4 m/hour at 1000 mD, 5 bar driving pressure).

### Known limitations

- The reactive transport coupling (deliverable 7 in the workflow) is
  not yet included; this initial release covers flow physics only.
- The Henry-law adapter is intentionally left at PFLOTRAN's default
  (bulk-air `EOSGasHenry_air`); species-specific O₂/N₂/CO₂/Ar
  solubility is accessed via the wrapper for the chemistry coupling,
  not through PFLOTRAN's main mass balance.
- Tested against PFLOTRAN's recent development branch. Other PFLOTRAN
  versions may require minor adjustments to file paths or
  `auxvar` field names.
