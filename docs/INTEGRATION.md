# PFLOTRAN integration guide

This document describes how to integrate `Air_EOS_module` and its
four supporting modules into a working PFLOTRAN build. The integration
is more involved than dropping in source files because PFLOTRAN's
GENERAL mode uses a hard-coded dispatcher pattern that calls CO₂-EOS
routines by name — we need to either:

  (a) replace those calls with calls into our wrapper when an "air"
      fluid type is selected, or
  (b) shadow the CO₂ EOS in a separate compilation unit and select
      between them at the input-deck level.

Approach (a) is what this guide describes. It is the cleaner, longer-
lived solution, but requires modest patches in three PFLOTRAN source
files.

Start with this guide as a *checklist*, not a literal patch — the
line numbers below reflect a PFLOTRAN source tree consistent with my
training-data snapshot. The exact line numbers in the live
`master`/`dev` branch will differ; the logical structure should be
nearly identical.

## Step 0. Sanity check the standalone build

Before touching PFLOTRAN, confirm the standalone tests pass on your
machine:

```bash
cd ~/air_eos_pr
make clean && make all
make test sweep henry watervapour properties wrapper
make crosscheck    # optional Python independent verification
```

You should see all six test executables produce uniform "PASS"
output. If anything fails here, fix it before touching PFLOTRAN.

## Step 1. Copy source files into PFLOTRAN

From PFLOTRAN's top-level directory:

```bash
PFLOTRAN_AIR_EOS=$HOME/air_eos_pr
cp $PFLOTRAN_AIR_EOS/air_eos_pr_module.F90       src/pflotran/
cp $PFLOTRAN_AIR_EOS/air_henry_module.F90        src/pflotran/
cp $PFLOTRAN_AIR_EOS/air_watervapour_module.F90  src/pflotran/
cp $PFLOTRAN_AIR_EOS/air_properties_module.F90   src/pflotran/
cp $PFLOTRAN_AIR_EOS/air_eos_module.F90          src/pflotran/
```

Five new files. Each contains the preprocessor block

```fortran
#ifndef PFLOTRAN_INTEGRATION
#define PetscReal      real(kind=8)
...
#endif
```

so all you need is to ensure `-DPFLOTRAN_INTEGRATION` reaches the
compiler. PFLOTRAN typically adds `-cpp` and `-D` flags via PETSc
configuration; check `makefile` and `petscvariables` to confirm.

## Step 2. Add files to the build

Find PFLOTRAN's source list, typically near the top of
`src/pflotran/pflotran_object_files.txt` or in the makefile. Add the
five new objects in dependency order:

```
air_eos_pr_module.o
air_henry_module.o
air_watervapour_module.o
air_properties_module.o
air_eos_module.o
```

Order matters because `air_eos_module.F90` `use`s the other four. If
the build uses automatic dependency tracking via `makedepf90`, you
can ignore order — the dependencies are detected from `use`
statements.

## Step 3. Register the new fluid type in the deck parser

Find PFLOTRAN's general-mode option parser. In recent PFLOTRAN
versions this lives in `src/pflotran/general.F90` around the
`GeneralReadOptions` subroutine. Find the case statement that
handles existing fluid types:

```fortran
case('AIR_WATER')
   option%nphase = 2
   option%nflowdof = 3
   ...
case('CO2_WATER')
   ...
```

The `AIR_WATER` keyword probably already exists with a simple ideal-
gas treatment. You have two clean choices:

**Option 3a.** Override `AIR_WATER` to use our wrapper. Add a new
internal flag like `option%use_high_pressure_air_eos = .true.` set
when entering this case. Then anywhere `AIR_WATER` mode currently
uses ideal-gas air properties (search for `IDEAL_GAS` in
`general_aux.F90`), guard with that flag and dispatch to
`Air_EOS_module` instead.

**Option 3b.** Add a parallel keyword `AIR_WATER_HP` (or similar)
for "high-pressure air":

```fortran
case('AIR_WATER_HP')
   option%nphase = 2
   option%nflowdof = 3
   option%use_high_pressure_air_eos = .true.
```

Option 3b is safer because it leaves existing deck behaviour
unchanged. Recommend 3b for an integration sandbox; switch to 3a
once everything works.

## Step 4. Patch the gas-phase property dispatcher

Open `src/pflotran/general_aux.F90`. Find the routine that updates
gas-phase auxiliary properties — typically called
`GeneralAuxUpdateGasPhase` or similar. Inside that routine you will
find calls to compute gas density, viscosity, enthalpy, and fugacity
from the current pressure-temperature state.

Replace those calls (gated on the high-pressure-air flag) with a
single call into our wrapper:

```fortran
#include "petsc/finclude/petscsys.h"
use Air_EOS_module

! ... existing local variables ...

PetscReal :: rho_gas, mu_gas, fug_coeff, h_dep
PetscReal :: drho_dP, drho_dT, dmu_dP, dmu_dT
PetscReal :: dfug_dP, dfug_dT, dhdep_dP, dhdep_dT
PetscErrorCode :: ierr

! ...

if (option%use_high_pressure_air_eos) then
  call AirEOSGetGasProperties(auxvar%pres(GAS_PHASE),                  &
                              auxvar%temp,                             &
                              rho_gas, mu_gas, fug_coeff, h_dep,       &
                              drho_dP, drho_dT,                        &
                              dmu_dP, dmu_dT,                          &
                              dfug_dP, dfug_dT,                        &
                              dhdep_dP, dhdep_dT, ierr)
  if (ierr /= 0) then
    option%status = 'AirEOS failure'
    return
  endif
  ! Store into auxvar in the same fields the existing ideal-gas
  ! path uses:
  auxvar%den(GAS_PHASE)     = rho_gas / AIR_MW_KG  ! convert kg/m^3 -> mol/m^3 if PFLOTRAN uses molar
  auxvar%den_kg(GAS_PHASE)  = rho_gas
  auxvar%visc(GAS_PHASE)    = mu_gas
  ! ... etc. for all fields PFLOTRAN GENERAL mode tracks ...
else
  ! existing ideal-gas treatment
endif
```

The exact field names depend on the PFLOTRAN version. Look at how
the existing CO₂ EOS integration writes into `auxvar` — copy that
pattern.

## Step 5. Patch the aqueous-phase composition dispatcher

Same `general_aux.F90`, different routine — typically
`GeneralAuxUpdateLiquidPhase` or the routine that computes dissolved
gas mole fractions. Replace the existing Henry's-law call with:

```fortran
if (option%use_high_pressure_air_eos) then
  call AirEOSGetAirSolubility(auxvar%pres(GAS_PHASE),                   &
                              auxvar%temp,                              &
                              auxvar%salinity,                          &
                              y_gas_air_composition,                    &
                              x_aq_4species,                            &
                              dxaq_dP, dxaq_dT, dxaq_dm, ierr)
  ! x_aq_4species(HENRY_O2) is the O2 aqueous mole fraction
  ! similarly for N2, CO2, Ar
endif
```

The composition vector `y_gas_air_composition` should be:

```fortran
y_gas_air_composition(HENRY_O2)  = 0.20946d0
y_gas_air_composition(HENRY_N2)  = 0.78084d0
y_gas_air_composition(HENRY_CO2) = 4.21d-4
y_gas_air_composition(HENRY_AR)  = 0.00934d0
```

unless your deck specifies a non-air composition.

## Step 6. Verify the integration

PFLOTRAN regression-test convention: drop input decks in
`regression_tests/general/air_caes/`. Start with the simplest
possible case — single cell, no chemistry, just gas in a pressurized
reservoir — and verify that:

1. The simulation runs without errors.
2. The reported gas density matches the standalone module output
   when called at the same (P, T).
3. The reported viscosity matches.
4. The reported Henry's-law dissolved O₂ matches.

If 1–4 all pass for a single cell, escalate to a 1D column with one
injection-production cycle.

Deliverable 6 will provide the actual input deck templates. For now,
the standalone test executables can produce reference values:

```bash
./test_air_eos_module    # prints reference values at 80 bar, 40 C
```

## Step 7. Call AirEOSVerify at startup

Optional but recommended: in PFLOTRAN's main driver, after
input-deck parsing but before the time loop, call:

```fortran
if (option%use_high_pressure_air_eos) then
  call AirEOSVerify(.true., ierr)
  if (ierr /= 0) then
    write(*,*) 'Air EOS self-check failed, ierr =', ierr
    call PetscFinalize(petsc_ierr)
    stop
  endif
endif
```

This catches any compile/link issues where the wrapper was
incorporated but is mis-linked.

## Common pitfalls

1. **Wrong gas constant.** PFLOTRAN uses an `IDEAL_GAS_CONSTANT`
   defined in `PFLOTRAN_Constants_module`. Our modules use it when
   compiled with `-DPFLOTRAN_INTEGRATION`. If you see a small
   ~ 1e-6 discrepancy between standalone and PFLOTRAN-integrated
   results, this is the cause — both should match to bit precision
   if the same constant is used.

2. **Mass fraction vs. mole fraction conventions.** PFLOTRAN's
   `auxvar%xmol_*` is mole fraction; `auxvar%xmass_*` is mass
   fraction. The wrapper returns mole fractions. Convert via
   `Air_MW = 28.96 g/mol` and `H2O_MW = 18.015 g/mol` as needed.

3. **Pressure conventions.** PFLOTRAN sometimes uses
   "Pa above atmospheric" and sometimes absolute Pa. The wrapper
   uses absolute Pa. If you get a numerical NaN near `P -> 0`,
   suspect a sign or offset confusion.

4. **Primary-variable switching during two-phase transitions.**
   PFLOTRAN's GENERAL mode switches primary variables when a phase
   appears or disappears. For CAES at moderate T (< 80 C), air is
   always in the supercritical/gas branch and water always has a
   liquid phase — the transitions are simpler than for CO₂-water.
   You may be able to disable some of the more complex switching
   logic and use a simpler liquid-saturated assumption. Test this
   carefully on a small problem first.

5. **PFLOTRAN's CO2_WATER vs your AIR_WATER_HP keyword.** If you
   accidentally pick `CO2_WATER` mode while integrating, the CO2
   Span-Wagner dispatcher will still be active. Always confirm
   `option%fluid_name == 'AIR_WATER_HP'` before claiming any
   integration result.

## When this gets stuck

Two strong external resources:

- **Glenn Hammond** (Sandia, PFLOTRAN lead): worth emailing before
  you start. He may already have an air-EOS branch sitting around
  unmerged. PFLOTRAN's GitHub issues are also a reasonable place to
  surface the integration intent.

- **PFLOTRAN-Dev mailing list** (https://groups.google.com/g/pflotran-dev):
  active, responsive to integration questions.

In either case, lead with what you've already done (verified
standalone EOS modules with full analytic derivatives) — that
position makes for a much shorter conversation than starting from
"how do I add a new EOS to PFLOTRAN".
