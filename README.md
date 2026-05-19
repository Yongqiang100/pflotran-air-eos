# pflotran-air-eos

High-pressure air equation-of-state modules for PFLOTRAN, built for
compressed-air energy storage (CAES) hydrogeochemistry simulations.

This repository contains:

- A set of standalone Fortran modules implementing a verified
  Peng-Robinson + Peneloux air EOS, Lemmon-Jacobsen viscosity,
  Henry's law solubility for O₂/N₂/CO₂/Ar, IAPWS-IF97 water-vapour
  partitioning, and NIST Shomate ideal-gas enthalpy
- Integration code that hooks these modules into PFLOTRAN's GENERAL
  mode via the function-pointer dispatcher in `eos_gas.F90`
- Regression test decks at increasing physical complexity
- A complete workflow document describing how the integration was
  performed (reproducible for other EOS additions)

## Status

| Component | Status |
|---|---|
| Standalone EOS modules | Verified against NIST/literature anchors |
| Python cross-validation | Bit-precision agreement with Fortran |
| PFLOTRAN integration | Successfully built and numerically validated |
| Single-cell verification | ρ_g matches standalone to 0.03 % |
| 1D flow verification | Newton converges in ~2 iter/step; gas front advances physically |
| Geochemistry coupling | In development (chemistry network design phase) |

## Motivation

Compressed-air energy storage (CAES) in saline aquifers introduces O₂
into anoxic groundwater containing reduced minerals (pyrite,
arsenopyrite, organic matter). The resulting oxidation reactions
mobilise trace metals — particularly arsenic — and can degrade
groundwater quality in adjacent fresh-water aquifers if cap-rock seals
leak.

Prior CAES reservoir modeling (Yang et al. 2024, Wang & Bauer 2019)
focuses on flow physics with simplified or absent chemistry, citing
minimal permeability/porosity feedback from reactions. This is correct
for reservoir engineering but does not answer the groundwater-quality
question: even with negligible feedback on flow, mobilized metals
can leak across cap-rock at concentrations harmful to drinking-water
quality.

To address this gap rigorously, reactive transport simulations need:

1. **An accurate high-pressure air EOS** (ideal gas is ~5–10 % off at
   CAES pressures, distorting fugacity-corrected Henry's law and
   therefore O₂ availability for oxidation)
2. **Henry's law for individual species** (O₂, N₂, Ar, CO₂) rather
   than bulk-air solubility, since the chemistry network operates on
   individual dissolved gases
3. **Water vapour partitioning** under the elevated total pressure
4. **Viscosity reflecting real gas density**, since the Darcy flux
   denominator changes ~10 % between atmospheric and CAES conditions

This repository provides those pieces, integrated into PFLOTRAN.

## Quick start

### Prerequisites

- Linux machine (Ubuntu, Debian, RHEL, etc.)
- `gfortran` 11+ with `-cpp` support
- (For PFLOTRAN integration) PETSc and PFLOTRAN built from source
  on the same machine

### Build and verify the standalone modules

```bash
git clone https://github.com/Yongqiang100/pflotran-air-eos.git
cd pflotran-air-eos/src
make
make test wrapper henry watervapour properties
```

You should see uniform `PASS` output across all six test executables.
Each test runs in well under a second.

### Integrate into PFLOTRAN

Detailed steps in [`docs/WORKFLOW.md`](docs/WORKFLOW.md). Briefly:

1. Build PETSc (see PETSc docs)
2. Build PFLOTRAN baseline
3. Copy the five `.F90` modules from `src/` into PFLOTRAN's
   `src/pflotran/`
4. Add entries to `pflotran_object_files.txt` (under `eos_obj`)
   and `pflotran_dependencies.txt`
5. Apply the patches documented in
   [`pflotran_integration/`](pflotran_integration/) to add adapter
   routines and the `PR_PENELOUX_AIR` keyword
6. Rebuild PFLOTRAN
7. Run regression tests from [`regression_tests/`](regression_tests/)

### Activate the EOS in an input deck

```
SUBSURFACE

EOS GAS
  DENSITY PR_PENELOUX_AIR
END

# ... rest of standard PFLOTRAN GENERAL-mode deck ...

END_SUBSURFACE
```

## Verified accuracy

At a representative CAES state (80 bar, 40 °C, 0.5 mol/kg NaCl brine):

| Property | Reference value | Source |
|---|---|---|
| Gas density | 87.95 kg/m³ | This work, vs NIST sweep |
| Gas viscosity | 2.06 × 10⁻⁵ Pa·s | Lemmon-Jacobsen 2004 |
| Fugacity coefficient | 1.003 | This work |
| Water vapour fraction | 9.6 × 10⁻⁴ | IAPWS-IF97 + Poynting |
| Dissolved O₂ | 14.6 mol/m³ (466 mg/L) | Sander 2015 + Setschenow |

PFLOTRAN-integrated values match the standalone reference to within
0.03 % for density and within numerical roundoff for viscosity.

## Repository structure

```
pflotran-air-eos/
├── README.md                       This file
├── LICENSE                         MIT
├── CITATION.cff                    Academic citation metadata
├── CHANGELOG.md                    Version history
├── docs/
│   ├── WORKFLOW.md                 Complete integration workflow
│   ├── INTEGRATION.md              PFLOTRAN integration guide
│   └── methods.md                  Publishable methods description
├── src/                            Standalone Fortran modules
│   ├── Makefile
│   ├── README.md
│   ├── air_eos_pr_module.F90       Peng-Robinson + Peneloux EOS
│   ├── air_henry_module.F90        Henry's law for O₂/N₂/CO₂/Ar
│   ├── air_watervapour_module.F90  IAPWS-IF97 + Poynting + activity
│   ├── air_properties_module.F90   Viscosity (L-J) + ideal-gas h
│   ├── air_eos_module.F90          PFLOTRAN-facing wrapper
│   ├── test_*.F90                  Test drivers (one per module)
│   └── cross_validate.py           Independent Python verification
├── pflotran_integration/
│   ├── README.md
│   └── integration_patch.txt       Source-code patches for eos_gas.F90
│                                   and eos.F90
└── regression_tests/
    ├── README.md
    ├── level1_single_cell.in       Single-cell EOS dispatch check
    ├── level2_pressure_relax.in    Single-cell pressure transient
    ├── level3_1d_injection.in      1D column gas injection
    └── level5_2d_yang_singlecycle.in   2D radial (Yang et al. geometry)
```

## Documentation

- [`docs/WORKFLOW.md`](docs/WORKFLOW.md) — Complete recipe for
  reproducing or extending the integration. Includes lessons learned
  and common pitfalls.
- [`docs/INTEGRATION.md`](docs/INTEGRATION.md) — Step-by-step PFLOTRAN
  integration guide.
- [`docs/methods.md`](docs/methods.md) — Publication-style description
  of the formulations, suitable for citing in a paper.

## Scientific references

The physical formulations used here:

- Peng, D.Y., Robinson, D.B. (1976). A new two-constant equation of
  state. *Ind. Eng. Chem. Fund.* 15, 59–64.
- Péneloux, A., Rauzy, E., Fréze, R. (1982). A consistent correction
  for Redlich-Kwong-Soave volumes. *Fluid Phase Equilib.* 8, 7–23.
- Lemmon, E.W., Jacobsen, R.T., Penoncello, S.G., Friend, D.G. (2000).
  Thermodynamic properties of air. *J. Phys. Chem. Ref. Data* 29,
  331–385.
- Lemmon, E.W., Jacobsen, R.T. (2004). Viscosity and thermal
  conductivity equations for nitrogen, oxygen, argon, and air. *Int.
  J. Thermophysics* 25, 21–69.
- Wagner, W., Pruß, A. (2002). The IAPWS formulation 1995 for the
  thermodynamic properties of ordinary water substance for general
  and scientific use. *J. Phys. Chem. Ref. Data* 31, 387–535.
- Sander, R. (2015). Compilation of Henry's law constants. *Atmos.
  Chem. Phys.* 15, 4399–4981.
- Weisenberger, S., Schumpe, A. (1996). Estimation of gas solubilities
  in salt solutions at temperatures from 273 K to 363 K. *AIChE J.*
  42, 298–300.

Prior CAES modeling literature this work builds on:

- Yang, Z., Lu, C., Liu, B., Yu, Y. (2024). Multi-cycle simulation of
  compressed air energy storage in aquifers. *J. Energy Storage* 99,
  113202.
- Wang, B., Bauer, S. (2019). Pressure response of large-scale
  compressed air energy storage in porous formations. *Appl.
  Geochem.* 102, 171–185.
- Medeiros, M. et al. (2018). [PG&E pilot field study, EPRI report].

## Citation

If you use this code in academic work, please cite:

```bibtex
@software{pflotran_air_eos,
  author       = {Yongqiang Chen},
  title        = {pflotran-air-eos: High-pressure air EOS for PFLOTRAN
                  in compressed-air energy storage applications},
  year         = {2026},
  url          = {https://github.com/Yongqiang100/pflotran-air-eos},
  version      = {0.1.0}
}
```

A journal publication describing the chemistry application is in
preparation.

## Contributing

Issues and pull requests welcome. For substantive changes (new EOS
variants, additional gas mixtures, alternative mixing rules), please
open an issue first to discuss the design.

## License

MIT — see [`LICENSE`](LICENSE). You may use this code freely in
academic or commercial work, with attribution.

## Acknowledgments

PFLOTRAN's GENERAL-mode architecture (the function-pointer
dispatcher in `eos_gas.F90`) made this integration considerably
cleaner than alternative approaches would have been. Thanks to
the PFLOTRAN development team for that design.

The reaction-network design for the chemistry coupling draws on
PHREEQC's Dzombak-Morel surface-complexation database.
