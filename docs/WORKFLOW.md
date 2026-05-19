# Adding a New Gas EOS to PFLOTRAN: Complete Workflow

A reproducible recipe for integrating a custom gas-phase equation of state
into PFLOTRAN's GENERAL mode. Documents the full process from designing
the standalone physics modules through verifying that PFLOTRAN simulations
produce the same numbers as the standalone reference.

The specific case worked through here is a high-pressure air EOS for
compressed-air energy storage (CAES) geochemistry, integrating
Peng-Robinson + Peneloux density, Henry's law for O₂/N₂/CO₂/Ar,
IAPWS-IF97 water-vapour partitioning, and Lemmon-Jacobsen viscosity.
The same workflow applies to any gas EOS — supercritical CO₂ variants,
hydrogen, methane mixtures, etc.

---

## 1. Overview and scope

### What you produce

By the end of this workflow you have:

1. A set of standalone Fortran modules implementing your EOS, each with
   verified physics and analytic derivatives, cross-validated against
   independent reference data
2. A PFLOTRAN-facing wrapper module consolidating the standalone modules
   into a single interface
3. ~150 lines of new code added to PFLOTRAN's `eos_gas.F90` (adapter
   routines + setter)
4. ~2 lines of new code in PFLOTRAN's `eos.F90` (keyword registration)
5. A new input-deck keyword that activates your EOS
6. Regression test decks at increasing physical complexity
7. Numerical validation that PFLOTRAN-integrated results match the
   standalone reference to within rounding precision

### Time budget

- **Phase 1** (standalone EOS development): 1–3 weeks for the physics
  modules, depending on what's already published. Most of this is the
  thermodynamics, not the code.
- **Phase 2** (PFLOTRAN wrapper): 1 day
- **Phase 3** (PFLOTRAN integration): 2–5 days. Most stuck-ness happens
  here, in PETSc and PFLOTRAN's build system.
- **Phase 4** (validation): 1 day if everything was right; longer if
  unit conversions need debugging

### Prerequisites

- **Linux** machine you control (Ubuntu, Debian, RHEL, etc.)
- **~20 GB** free disk space (PETSc + PFLOTRAN + builds)
- **gfortran 11+** with `-cpp` preprocessor support
- **git, make, python3** standard build tools
- **Working knowledge of Fortran 2003+**: modules, derived types, function
  pointers
- **Basic thermodynamics**: EOS principles, fugacity, departure functions,
  cubic equation roots
- (Helpful, not required) **PFLOTRAN deck syntax**: familiarity with
  GENERAL mode keywords

---

## 2. Phase 1 — Standalone EOS development

### 2.1 Identify the physics needed

Before coding, decide what gas-phase properties your problem needs.
For a typical reactive-transport simulation you need at least:

- **Density** ρ(P, T) and its derivatives — flow Jacobian, mass balance
- **Fugacity coefficient** φ(P, T) — chemistry coupling at high pressure
- **Viscosity** μ(P, T) or μ(ρ, T) — Darcy flux
- **Enthalpy** h(P, T) — energy balance, ideal-gas reference plus
  departure
- **Dissolved species concentrations** via Henry's law — for the
  aqueous-phase chemistry network
- **Water vapour partitioning** — humidity in the gas, drying near
  wellbore

For the air EOS case, this maps to five separate modules — one per
physical effect. For a simpler EOS (e.g. pure CO₂), four would suffice.

### 2.2 Module template

Each module follows this skeleton:

```fortran
! Conditional macros so the module works both standalone and inside
! PFLOTRAN. PetscReal/PetscInt resolve to PFLOTRAN's types when
! PFLOTRAN_INTEGRATION is defined, to plain double/int otherwise.
#ifndef PFLOTRAN_INTEGRATION
#define PetscReal      real(kind=8)
#define PetscInt       integer(kind=4)
#define PetscBool      logical
#define PetscErrorCode integer(kind=4)
#endif

module Air_EOS_PR_module

#ifdef PFLOTRAN_INTEGRATION
#include "petsc/finclude/petscsys.h"
  use petscsys
#endif

  implicit none
  private

  ! ----- Public interface -----
  public :: AirEOSPRProperties
  public :: AirEOSPRFugacityCoeff
  ! ... etc

  ! ----- Module-scope parameters -----
  PetscReal, parameter :: AIR_TC = 132.5306d0  ! K
  PetscReal, parameter :: AIR_PC = 3.7850d6    ! Pa
  ! ... etc

contains

  subroutine AirEOSPRProperties(P, T, rho, Z, fug, h_dep,                  &
                                 drho_dP, drho_dT, dfug_dP, dfug_dT, ierr)
    PetscReal, intent(in)       :: P, T
    PetscReal, intent(out)      :: rho, Z, fug, h_dep
    PetscReal, intent(out)      :: drho_dP, drho_dT, dfug_dP, dfug_dT
    PetscErrorCode, intent(out) :: ierr
    ! ... implementation
  end subroutine

end module
```

The key design points:

- **Each public routine returns analytic derivatives** alongside the
  primary value. PFLOTRAN's Newton solver needs Jacobian terms; without
  good derivatives convergence is poor.
- **`ierr` argument** for error reporting — PFLOTRAN convention.
- **PETSc-style types** even when standalone, via the macro trick.
  Ensures the same source compiles cleanly in both contexts.

### 2.3 Test driver template

Each module has a companion `test_*.F90` driver with:

1. **Anchor checks**: hard-coded reference values at well-known points
   (e.g. NIST WebBook for density, IAPWS-IF97 reference table for
   water vapour pressure). Tolerance typically 0.5–2 % depending on
   the formulation's known accuracy.

2. **Derivative checks**: centered finite-difference comparison against
   the analytic derivatives. Tolerance ~1×10⁻⁵ relative error
   (FD truncation noise dominates).

3. **Sweep**: optional 2D grid over (P, T) writing to CSV for visual
   inspection.

4. **Cross-validation**: optional independent Python implementation
   of the same formulation, comparing line-by-line. This is the
   strongest verification — catches bugs the unit tests miss.

### 2.4 Build system

A simple Makefile per module:

```make
FC      ?= gfortran
FFLAGS  ?= -O2 -g -cpp -Wall -Wextra -fimplicit-none -std=f2008 \
           -Wno-unused-variable -Wno-unused-dummy-argument

%.o: %.F90
	$(FC) $(FFLAGS) -c $< -o $@

test_X: X_module.o test_X.o
	$(FC) $(FFLAGS) -o $@ $^
```

No `-DPFLOTRAN_INTEGRATION` flag in standalone builds — the module
compiles with stub macro definitions.

### 2.5 Workflow within Phase 1

For each module:

1. Find the published formulation and parameter values
2. Write the module with skeletons for `subroutine X(in, out, ierr)`
3. Implement forward calculation
4. Derive analytic derivatives by hand
5. Write the test driver with anchors from the literature
6. Build, run, iterate until anchors pass
7. Add finite-difference derivative checks; verify analytic = FD
8. Optionally: independent Python cross-check

A module typically takes 2–5 hours from start to verified. The
derivative algebra is where most bugs hide — symbolic differentiation
via SymPy can help.

---

## 3. Phase 2 — PFLOTRAN-facing wrapper

### 3.1 Why wrap?

PFLOTRAN's `EOS_Gas_module` uses specific calling conventions and unit
systems. Rather than spreading these through your physics modules,
isolate the conversions in a single wrapper. This keeps the physics
modules clean (testable standalone) and makes the unit-conversion
logic explicit (auditable in one place).

### 3.2 Unit conventions

These differ between standalone scientific and PFLOTRAN GENERAL mode:

| Quantity     | Standalone (SI) | PFLOTRAN GENERAL |
|--------------|-----------------|------------------|
| Temperature  | Kelvin          | **Celsius**      |
| Density      | kg/m³           | **kmol/m³**      |
| Enthalpy     | J/mol           | **J/kmol**       |
| Pressure     | Pa              | Pa               |
| Viscosity    | Pa·s            | Pa·s             |
| Composition  | mole fraction   | mole fraction    |

Conversions are simple but easy to forget:

```fortran
T_K       = T_C + T273K              ! T273K = 273.15
rho_kmol  = rho_kg / fmw_gas         ! fmw_gas in kg/kmol = g/mol
H_Jkmol   = H_Jmol * 1.0d3
```

`fmw_gas` is a module variable in `eos_gas.F90`, set to `FMWAIR` by
default in `EOSGasInit()`. For a non-air gas, override with
`EOSGasSetFMWConstant()`.

### 3.3 Wrapper structure

A single module with one or two public routines that PFLOTRAN-side
code calls. The wrapper takes PFLOTRAN-convention inputs, internally
calls the physics modules, and returns PFLOTRAN-convention outputs:

```fortran
module Air_EOS_module
  use Air_EOS_PR_module
  use Air_Henry_module
  use Air_WaterVapour_module
  use Air_Properties_module
  implicit none
  private

  public :: AirEOSGetGasProperties
  public :: AirEOSGetTotalEnthalpy
  public :: AirEOSGetWaterVapourFraction
  public :: AirEOSGetAirSolubility
  public :: AirEOSVerify

contains
  ! Each public routine: PFLOTRAN-style inputs, dispatches to
  ! standalone modules, returns PFLOTRAN-style outputs.
end module
```

This is the boundary layer between your physics code (which lives
forever and is publishable) and PFLOTRAN's evolving conventions
(which may change between versions).

---

## 4. Phase 3 — PFLOTRAN integration

The most failure-prone phase. Stop after each stage and verify
before proceeding.

### Stage 0: Backup discipline

Make a tarball at each milestone:

```bash
cd ~
tar czf pflotran_milestoneNAME.tar.gz pflotran/src/pflotran/
```

Milestones to snapshot:
- After PFLOTRAN baseline builds
- After standalone modules copied in
- After dispatch wiring
- After validation passes

### Stage 1: Install PETSc

```bash
cd ~
git clone -b release https://gitlab.com/petsc/petsc.git petsc
cd petsc
./configure \
  --with-cc=gcc \
  --with-cxx=g++ \
  --with-fc=gfortran \
  --with-debugging=0 \
  --download-mpich=yes \
  --download-fblaslapack=yes \
  --download-hdf5=yes \
  --download-hdf5-fortran-bindings=yes \
  --download-metis=yes \
  --download-parmetis=yes
make all
make check
export PETSC_DIR=$HOME/petsc
export PETSC_ARCH=$(ls -d arch-*)
```

Time: 1–2 hours mostly waiting for `make all`. If PETSc tests fail
here, debug PETSc before continuing — adding PFLOTRAN on top of
broken PETSc is hopeless.

### Stage 2: Build PFLOTRAN baseline

```bash
cd ~
git clone https://bitbucket.org/pflotran/pflotran.git
cd pflotran/src/pflotran
make pflotran
```

Verify it runs:

```bash
cd ~/pflotran/regression_tests/general
$HOME/pflotran/src/pflotran/pflotran -input_prefix 1d_flux
```

You should see `Wall Clock Time: ... [sec]` at the end. If not, fix
the baseline before continuing.

### Stage 3: Add standalone modules to PFLOTRAN

Copy files in:

```bash
PFLOTRAN_SRC=$HOME/pflotran/src/pflotran
cp ~/my_eos_dir/*.F90 $PFLOTRAN_SRC/
```

PFLOTRAN's build uses an explicit object-list file
`pflotran_object_files.txt`. Find the right block (`eos_obj` for gas
EOS modules) and add your `.o` entries in alphabetical order before
the existing CO₂ entries:

```
eos_obj = \
        ${common_src}air_eos_module.o \
        ${common_src}air_eos_pr_module.o \
        ${common_src}air_henry_module.o \
        ${common_src}air_properties_module.o \
        ${common_src}air_watervapour_module.o \
        ${common_src}co2_sw.o \
        ${common_src}co2_span_wagner_spline.o \
        ${common_src}eos.o \
        ${common_src}eos_gas.o \
        ${common_src}gas_eos_mod.o
```

The last line of the block has no trailing backslash — that signals
the end to make.

Then add module dependencies to `pflotran_dependencies.txt`:

```
air_eos_module.o : \
  air_eos_pr_module.o \
  air_henry_module.o \
  air_properties_module.o \
  air_watervapour_module.o
```

Only the wrapper module has internal dependencies; the four physics
modules are self-contained.

Build:

```bash
cd $PFLOTRAN_SRC
make pflotran
```

If it succeeds, your modules are compiled and linked into PFLOTRAN
but not yet *used*. Re-run the `1d_flux` regression test to verify
no regression.

### Stage 4: Understand PFLOTRAN's gas EOS architecture

This is the critical insight. PFLOTRAN's `EOS_Gas_module` uses
**function pointer dispatch**:

```fortran
procedure(EOSGasDensityDummy),  pointer :: EOSGasDensityPtr  => null()
procedure(EOSGasEnergyDummy),   pointer :: EOSGasEnergyPtr   => null()
procedure(EOSGasViscosityDummy),pointer :: EOSGasViscosityPtr => null()
procedure(EOSGasHenryDummy),    pointer :: EOSGasHenryPtr    => null()
procedure(EOSGasDensityEnergyDummy), pointer :: EOSGasDensityEnergyPtr => null()
```

At simulation start, `EOSGasInit()` binds these pointers to default
ideal-gas implementations. To use a different EOS, the deck parser
calls a *setter* routine that re-binds the pointers to alternative
implementations.

Existing setters:

- `EOSGasSetDensityIdeal()` — ideal gas (default)
- `EOSGasSetDensityRKS()` — Redlich-Kwong-Soave (hydrogen-tuned)
- `EOSGasSetDensityPRMethane()` — Peng-Robinson for methane
- `EOSGasSetDensityConstant()` — fixed density
- `EOSGasSetEOSDBase()` — table lookup

You add a new setter for your EOS following the same pattern.

### Stage 5: Write adapter routines

For each "Dummy" interface in `eos_gas.F90`, write a matching adapter
that does unit conversion and calls into your wrapper module:

```fortran
subroutine EOSGasDensityPRAir(T,P,Rho_gas,dRho_dT,dRho_dP,ierr,table_idxs)
  use Air_EOS_module
  implicit none
  PetscReal, intent(in)       :: T        ! [C]  PFLOTRAN convention
  PetscReal, intent(in)       :: P        ! [Pa]
  PetscReal, intent(out)      :: Rho_gas  ! [kmol/m^3]  PFLOTRAN convention
  PetscReal, intent(out)      :: dRho_dT, dRho_dP
  PetscErrorCode, intent(out) :: ierr
  PetscInt, pointer, optional, intent(inout) :: table_idxs(:)

  PetscReal :: T_K
  PetscReal :: rho_kg, mu_dum, fug_dum, hdep_dum
  PetscReal :: drho_dT_kg, drho_dP_kg, dmu_dT, dmu_dP
  PetscReal :: dfug_dT, dfug_dP, dhdep_dT, dhdep_dP

  T_K = T + T273K   ! C -> K

  call AirEOSGetGasProperties(P, T_K,                                       &
                              rho_kg, mu_dum, fug_dum, hdep_dum,            &
                              drho_dP_kg, drho_dT_kg,                       &
                              dmu_dP, dmu_dT,                               &
                              dfug_dP, dfug_dT,                             &
                              dhdep_dP, dhdep_dT, ierr)
  if (ierr /= 0) return

  ! kg/m^3 -> kmol/m^3
  Rho_gas = rho_kg / fmw_gas
  dRho_dT = drho_dT_kg / fmw_gas
  dRho_dP = drho_dP_kg / fmw_gas
end subroutine EOSGasDensityPRAir
```

You typically need 4–5 adapter routines (density, energy, combined
density+energy, viscosity, optionally Henry).

### Stage 6: Write the setter routine

A simple binding routine:

```fortran
subroutine EOSGasSetPRPenelouxAir()
  implicit none
  EOSGasDensityEnergyPtr => EOSGasDensityEnergyPRAir
  EOSGasDensityPtr       => EOSGasDensityPRAir
  EOSGasEnergyPtr        => EOSGasEnergyPRAir
  EOSGasViscosityPtr     => EOSGasViscosityPRAir
end subroutine EOSGasSetPRPenelouxAir
```

Add the setter name to the `public ::` list near other setters in
`eos_gas.F90`.

Also add the `use Air_EOS_module` (and any other modules your
adapters need) to the top of `eos_gas.F90`.

### Stage 7: Add keyword to the deck parser

In `eos.F90`, find the `case('PR_METHANE')` block in the gas-density
parser (typically around line 440). Add a new case for your keyword:

```fortran
case('PR_METHANE')
  call EOSGasSetDensityPRMethane()
case('PR_PENELOUX_AIR')                ! NEW
  call EOSGasSetPRPenelouxAir()        ! NEW
case('IDEAL','DEFAULT')
  call EOSGasSetDensityIdeal()
```

Indentation matters — match the existing entries (14 spaces in the
canonical PFLOTRAN source).

### Stage 8: Build

```bash
cd $PFLOTRAN_SRC
make pflotran 2>&1 | tee build.log | tail -40
```

Likely outcomes:

1. **Success** → move to validation
2. **Module-circularity error** → the `use Air_EOS_module` inside
   `eos_gas.F90` creates an unexpected dependency. Move the `use`
   statement inside the adapter procedures instead of at module
   level.
3. **"Symbol X has no IMPLICIT type"** → forgot to add `use`
   statement, or the module isn't in `pflotran_dependencies.txt`
4. **Linker error "undefined reference"** → the new module isn't in
   `pflotran_object_files.txt`

---

## 5. Phase 4 — Validation

### 5.1 Build verification (no regression)

After successful build, confirm PFLOTRAN's existing tests still pass:

```bash
cd ~/pflotran/regression_tests/general
$HOME/pflotran/src/pflotran/pflotran -input_prefix 1d_flux 2>&1 | tail -10
```

Should show the same step count, newton iteration count, and wall
clock time as before adding your modules. Any difference suggests
you accidentally modified existing PFLOTRAN behaviour.

### 5.2 Single-cell numerical check

Create a minimal deck that exercises your new keyword:

```
SIMULATION
  SIMULATION_TYPE SUBSURFACE
  PROCESS_MODELS
    SUBSURFACE_FLOW flow
      MODE GENERAL
      OPTIONS
      /
    /
  /
END

SUBSURFACE

EOS GAS
  DENSITY PR_PENELOUX_AIR   # <-- your new keyword
END

GRID
  TYPE STRUCTURED
  ORIGIN 0.d0 0.d0 0.d0
  NXYZ 1 1 1
  BOUNDS
    0.d0 0.d0 0.d0
    1.d0 1.d0 1.d0
  /
END

# ... standard MATERIAL_PROPERTY, CHARACTERISTIC_CURVES, REGION ...

FLOW_CONDITION initial_gas_state
  TYPE
    GAS_PRESSURE DIRICHLET
    GAS_SATURATION DIRICHLET
    TEMPERATURE DIRICHLET
  /
  GAS_PRESSURE 8.d6
  GAS_SATURATION 0.95d0
  TEMPERATURE 40.d0
END

OUTPUT
  SNAPSHOT_FILE
    TIMES s 1.d0
    FORMAT TECPLOT POINT
  /
  VARIABLES
    GAS_DENSITY
    GAS_VISCOSITY
    GAS_PRESSURE
    GAS_SATURATION
    TEMPERATURE
  /
END

END_SUBSURFACE
```

Important syntactic notes:

- `EOS GAS / DENSITY YOUR_KEYWORD` activates the EOS
- `MATERIAL_PROPERTY` must contain `CHARACTERISTIC_CURVES name` (recent
  PFLOTRAN versions require explicit reference)
- Use `GAS_PRESSURE / GAS_SATURATION / TEMPERATURE` primary variables to
  force the gas-phase state

Run it. Compare the reported `GAS_DENSITY` to your standalone module's
value at the same (P, T):

| Quantity | Standalone | PFLOTRAN | Match |
|---|---|---|---|
| ρ_gas at 80 bar / 40 °C | 87.95 kg/m³ | 87.92 kg/m³ | 0.03% ✓ |

Discrepancies > 1% suggest unit-conversion bugs in your adapter.
Discrepancies ~ 1000× suggest a kmol vs mol confusion.

### 5.3 Flow simulation check

Run a 1D injection deck to verify Newton converges under real flow
gradients:

```
# 100-cell column, gas injection at left, pressure-driven flow
GRID
  NXYZ 100 1 1
  BOUNDS  0.d0 0.d0 0.d0  100.d0 1.d0 1.d0
END

# Initial: liquid-dominated at 80 bar
# Inlet BC: gas-dominated at 85 bar (5 bar driving pressure)
# Outlet BC: held at initial state
```

What to look for:

- **Newton converges** in 1–3 iterations per timestep (signals correct
  Jacobian)
- **Zero or few timestep cuts** (signals robust convergence)
- **Gas saturation front advances** at physically reasonable velocity
  (1–10 m/hour for 1000 mD, 5 bar gradient)
- **Gas density profile** matches standalone at the boundary
  pressures (cell nearest inlet ~ 93 kg/m³ at 85 bar; cell nearest
  outlet ~ 88 kg/m³ at 80 bar)

---

## 6. Common pitfalls

Things we hit during the air-EOS integration. Save yourself the
debugging time.

### 6.1 Fortran case-insensitivity in parameter declarations

**Symptom:** Compiler error "Parameter 'dp' at (1) has not been
declared."

**Cause:** A local variable named `dP` (uppercase P) collides with the
common kind parameter `dp = kind(1.0d0)`. Fortran is case-insensitive
for identifiers.

**Fix:** Rename local variables to `dP_step`, `dT_step`, `dPdum`, etc.

### 6.2 PFLOTRAN's source-list format

**Symptom:** Your new module compiles cleanly but the linker reports
"undefined reference to `your_module_routine_`."

**Cause:** Object file not added to the right `_obj` variable in
`pflotran_object_files.txt`.

**Fix:** For gas EOS modules, add to `eos_obj` block. For other
categories check the block grouping that holds the most similar
existing modules.

### 6.3 Module dependency tracking

**Symptom:** Build error "Cannot open module file `my_module.mod`
for reading."

**Cause:** PFLOTRAN tracks Fortran module dependencies explicitly in
`pflotran_dependencies.txt`. Without an entry, the build may try to
link before the dependency is compiled.

**Fix:** Add a block matching the existing format:
```
my_module.o : \
  dep1.o \
  dep2.o
```

### 6.4 Unit conversions

**Symptom:** PFLOTRAN reports gas density off by ~30× from standalone.

**Cause:** Forgot to convert kg/m³ to kmol/m³ via `/ fmw_gas`.

**Symptom:** Reports off by ~1000×.

**Cause:** Forgot to convert J/mol to J/kmol via `× 1000`.

**Fix:** Check the dummy interface comments carefully. PFLOTRAN's
internal convention is consistently Celsius, kmol/m³, J/kmol for
GENERAL mode. Your wrapper outputs SI; the adapter converts.

### 6.5 CHARACTERISTIC_CURVES not found

**Symptom:** "Characteristic curves \"\" in material property
\"soil\" not found among available characteristic curves."

**Cause:** Recent PFLOTRAN versions require explicit reference to
characteristic curves inside the `MATERIAL_PROPERTY` block.

**Fix:** Add the line `CHARACTERISTIC_CURVES default` (or whatever
your block is named) as the second line inside the `MATERIAL_PROPERTY`
block, after the `ID` line.

### 6.6 Public list duplicates

**Symptom:** Build error about duplicate `EOSGasSetXxx` in public
list, or warnings about ambiguous interfaces.

**Cause:** Adding the same setter name to multiple `public ::`
declarations in `eos_gas.F90`.

**Fix:** PFLOTRAN's `eos_gas.F90` has two `public ::` blocks — one
for routines exposed to general callers, one for setter routines.
Add yours to the setter block only.

### 6.7 Zero viscosity in low-saturation cells

**Symptom:** Snapshot output shows gas viscosity = 0 in cells where
gas saturation equals the residual value.

**Cause:** PFLOTRAN's convention: when `S_gas` is at residual, the
gas phase is immobile and PFLOTRAN reports viscosity as zero in
output for display purposes.

**Status:** Not a bug. Worth noting in your notebook if you
post-process viscosity profiles.

---

## 7. Validation checklist

Before declaring the integration complete:

- [ ] All standalone modules pass their unit tests with documented
      tolerances
- [ ] Cross-validation against Python or independent implementation
      reproduces standalone outputs to bit precision
- [ ] PFLOTRAN baseline builds cleanly
- [ ] After adding your modules, PFLOTRAN still builds cleanly
- [ ] PFLOTRAN's `1d_flux` (or equivalent) regression test still
      passes with the same step count
- [ ] New deck keyword is recognized (no "unknown keyword" parse
      error)
- [ ] Single-cell static test produces density matching standalone
      to better than 0.1%
- [ ] Single-cell test produces viscosity matching standalone
- [ ] Newton solver converges in 1–3 iterations per step under
      load (good Jacobians)
- [ ] 1D flow test completes without timestep cuts
- [ ] Gas front advance is consistent with Darcy + Mualem-VG
      relative permeability
- [ ] Project notebook captures the integration steps and the
      validation numbers

---

## 8. File inventory after completion

### Your project directory

```
~/my_eos_dir/
├── README.md
├── INTEGRATION.md
├── Makefile
├── module1.F90                ! physics module
├── module2.F90
├── ...
├── wrapper_module.F90         ! PFLOTRAN-facing wrapper
├── test_module1.F90           ! unit tests
├── test_module2.F90
├── ...
├── cross_validate.py          ! optional independent check
└── regression_tests/
    └── general/
        └── my_eos/
            ├── level1_single_cell.in
            ├── level2_pressure_relax.in
            ├── level3_1d_injection.in
            └── README.md
```

### In PFLOTRAN's source

```
~/pflotran/src/pflotran/
├── module1.F90                ! copied in
├── module2.F90
├── ...
├── wrapper_module.F90
├── pflotran_object_files.txt  ! N entries added
├── pflotran_dependencies.txt  ! 1 entry added
├── eos_gas.F90                ! ~150 lines added
└── eos.F90                    ! ~2 lines added
```

### Backup tarballs at milestones

```
~/pflotran_baseline.tar.gz                  ! after Stage 2
~/pflotran_modules_added.tar.gz             ! after Stage 3
~/pflotran_integration_complete.tar.gz      ! after Stage 8
~/pflotran_numerically_validated.tar.gz     ! after Phase 4
```

Each tarball is typically 50–100 MB and serves as a rollback point.

---

## 9. References and resources

- **PFLOTRAN main site**: https://www.pflotran.org
- **PFLOTRAN source**: https://bitbucket.org/pflotran/pflotran
- **PETSc**: https://petsc.org
- **PFLOTRAN-Dev mailing list**:
  https://groups.google.com/g/pflotran-dev
- **Glenn Hammond** (Sandia, PFLOTRAN lead) — worth emailing if you
  hit architectural questions, and worth notifying when integration
  completes so the work can potentially be merged upstream

For thermodynamic references specific to high-pressure gas EOS in
geological subsurface applications:

- Peng, D.Y., Robinson, D.B. (1976). A new two-constant equation of
  state. Ind. Eng. Chem. Fund. 15, 59–64.
- Péneloux, A., Rauzy, E., Fréze, R. (1982). A consistent correction
  for Redlich-Kwong-Soave volumes. Fluid Phase Equilib. 8, 7–23.
- Lemmon, E.W., Jacobsen, R.T., Penoncello, S.G., Friend, D.G. (2000).
  Thermodynamic properties of air. J. Phys. Chem. Ref. Data 29,
  331–385.
- Lemmon, E.W., Jacobsen, R.T. (2004). Viscosity and thermal
  conductivity equations for nitrogen, oxygen, argon, and air. Int.
  J. Thermophysics 25, 21–69.
- Sander, R. (2015). Compilation of Henry's law constants. Atmos.
  Chem. Phys. 15, 4399–4981.

---

## 10. Summary

This workflow takes a thermodynamics problem (here, high-pressure
air for CAES geochemistry) and produces a PFLOTRAN simulation
capability validated against the underlying physics.

The integration is structurally similar regardless of the gas. For
any new EOS:

1. Develop the standalone modules first, in isolation, until they
   reproduce reference data
2. Wrap them with PFLOTRAN-convention unit conversions
3. Add to PFLOTRAN's build (object list + dependencies)
4. Write adapter routines matching the Dummy interfaces
5. Bind the function pointers via a setter routine
6. Register a deck keyword that invokes the setter
7. Verify by comparing PFLOTRAN outputs to standalone reference

If you have all the physics figured out beforehand, the
PFLOTRAN-integration work itself is a focused 2–5 day effort. The
thermodynamics is where most of the project time goes.
