# Installation and Reproduction Guide

Complete end-to-end guide to reproduce the pflotran-air-eos work
from a fresh Linux machine. Every step is sequenced so that failure
at any point gives a clear diagnostic, not a cryptic compile error
two hours later.

This guide supersedes the original WORKFLOW.md, which mixed
"original development" with "reproduction" and didn't pin versions.

---

## What you'll have at the end

- A working PETSc v3.24.5 installation pinned to a specific tag
- A working PFLOTRAN baseline build verified against its own regression tests
- The five EOS Fortran modules integrated into PFLOTRAN
- Verified bit-precision match between standalone modules and the
  PFLOTRAN-integrated EOS at 80 bar / 40 °C
- A 1D injection simulation producing physically correct results

## Total time budget

- Active hands-on time: 2–3 hours
- Waiting for compilation: 1–3 hours (PETSc and PFLOTRAN builds)
- Total elapsed: ~4–6 hours, can be split across multiple sessions

## Required hardware

- Linux machine: Debian 12+, Ubuntu 22.04+, or compatible
- Windows machines: install WSL2 with Debian or Ubuntu first
- 20 GB free disk space (PETSc + PFLOTRAN builds are large)
- 4 GB RAM minimum
- Multi-core CPU helps (parallel `make -j`)

## Required software

These must be installed before you begin. The guide installs them
in Step 1, but check first:

```bash
gfortran --version    # need 11.0 or newer
gcc --version         # need 11.0 or newer
git --version         # any recent version
make --version        # GNU Make 4+
python3 --version     # 3.8+
```

If any are missing, install them with your package manager (Step 1
covers this).

---

## Step 1 — System dependencies

Time: 5 minutes.

### 1.1 Install build essentials

On Debian/Ubuntu:

```bash
sudo apt update
sudo apt install -y \
  build-essential \
  gfortran \
  git \
  make \
  cmake \
  pkg-config \
  python3 \
  python3-pip \
  wget \
  curl \
  ca-certificates
```

### 1.2 Verify

```bash
gfortran --version | head -1
gcc --version | head -1
git --version
make --version | head -1
python3 --version
df -h ~ | tail -1
```

You should see:

- gfortran 11 or higher
- gcc 11 or higher
- git 2.x
- GNU Make 4.x
- Python 3.8+
- At least 20 GB free on your home filesystem

If any of these check fails, fix it before proceeding. Do not try to
build PETSc without verifying the compiler.

---

## Step 2 — Build PETSc

Time: 30–60 minutes of compilation.

PETSc is a numerical solver library that PFLOTRAN depends on. The
critical insight: PFLOTRAN tracks a specific PETSc release tag.
Using `-b release` or `-b main` will give you whatever PETSc
released most recently, which may not match PFLOTRAN's expectations
and will produce cryptic Fortran binding errors.

**Always pin PETSc to the version PFLOTRAN's documentation
recommends.** Check the current recommendation at:

https://www.pflotran.org/documentation/user_guide/how_to/installation/linux.html

As of May 2026, that version is `v3.24.5`.

### 2.1 Clone and pin PETSc

```bash
cd ~
git clone https://gitlab.com/petsc/petsc.git petsc
cd petsc

# CRITICAL: pin to PFLOTRAN-supported version
git fetch --tags
git checkout v3.24.5

# Verify the checkout succeeded
git describe --tags
```

Expected output of `git describe --tags`:

```
v3.24.5
```

If you see anything else (a commit hash, "no tag", a different
version), STOP. The remaining steps depend on this being correct.

### 2.2 Configure PETSc

**Before configuring**, ensure PETSc environment variables aren't
pre-set to invalid values. PETSc rejects an empty `PETSC_ARCH` and
chooses its own default if these are unset. Setting them comes
later (§2.5), after the build.

```bash
unset PETSC_DIR PETSC_ARCH
echo "PETSC_DIR=[$PETSC_DIR]  PETSC_ARCH=[$PETSC_ARCH]"
```

Expected output:

```
PETSC_DIR=[]  PETSC_ARCH=[]
```

If either variable shows a value here, your shell still has them
exported (perhaps from a previous build attempt or from `.bashrc`).
The `unset` clears them for this shell session.

Now run PFLOTRAN's exact official configure command:

```bash
cd ~/petsc

./configure \
  --COPTFLAGS='-O3' \
  --CXXOPTFLAGS='-O3' \
  --FOPTFLAGS='-O3 -Wno-unused-function' \
  --with-debugging=no \
  --download-mpich=yes \
  --download-hdf5=yes \
  --download-hdf5-fortran-bindings=yes \
  --download-fblaslapack=yes \
  --download-metis=yes \
  --download-parmetis=yes
```

If configure fails with `PETSC_ARCH cannot be empty string`, you
didn't run the `unset` step above — go back and run it.

This will run for 5–15 minutes, downloading and configuring
MPICH, HDF5 (with Fortran bindings — this is essential), BLAS,
LAPACK, METIS, and ParMETIS. Watch for "Configure stage complete"
at the end.

If configure fails (rare on standard systems), the most common
causes are missing system packages (Step 1) or network issues.

### 2.3 Build PETSc

PETSc's configure prints the exact `make` line at the end. It
looks like:

```
xxx=========================================================================xxx
Configure stage complete. Now build PETSc libraries with:
make PETSC_DIR=/home/debian/petsc PETSC_ARCH=arch-linux-c-opt all
xxx=========================================================================xxx
```

Run that line (yours may differ slightly):

```bash
make PETSC_DIR=$HOME/petsc PETSC_ARCH=arch-linux-c-opt all
```

This takes 30–60 minutes depending on CPU. The output is verbose;
that's normal. Wait for it to finish without errors.

### 2.4 Verify PETSc

```bash
make PETSC_DIR=$HOME/petsc PETSC_ARCH=arch-linux-c-opt check
```

You should see:

```
Running check examples to verify correct installation
...
Completed test examples
```

If `make check` fails, do not proceed. Common causes:

- "Cannot find libpetsc.so": the build silently failed. Re-read
  the `make all` output for the actual error.
- "MPI tests failed": MPICH installation is broken. Try
  `--download-openmpi=yes` instead of mpich in 2.2.

### 2.5 Set environment variables

```bash
export PETSC_DIR=$HOME/petsc
export PETSC_ARCH=arch-linux-c-opt
echo "PETSC_DIR=$PETSC_DIR  PETSC_ARCH=$PETSC_ARCH"
```

Verify these are set. To make them persistent across shell
sessions, add to your `~/.bashrc`:

```bash
cat >> ~/.bashrc << 'BASHRC_EOF'

# PETSc for PFLOTRAN
export PETSC_DIR=$HOME/petsc
export PETSC_ARCH=arch-linux-c-opt
BASHRC_EOF
```

Open a new terminal and confirm:

```bash
echo "PETSC_DIR=$PETSC_DIR  PETSC_ARCH=$PETSC_ARCH"
```

Both should be set.

---

## Step 3 — Build PFLOTRAN baseline

Time: 30–60 minutes of compilation.

### 3.1 Clone PFLOTRAN

```bash
cd ~
git clone https://bitbucket.org/pflotran/pflotran
cd pflotran
git describe --tags --always
```

The `git describe` output tells you what commit/tag you're on.
Note this for your records.

### 3.2 Build PFLOTRAN

```bash
cd ~/pflotran/src/pflotran
make pflotran 2>&1 | tee /tmp/pflotran_baseline_build.log | tail -20
```

This takes 30–60 minutes. The output ends with linking — look for:

```
...
g++ ... -o pflotran ...
```

Success indicator: a `pflotran` executable in
`~/pflotran/src/pflotran/`.

```bash
ls -lh ~/pflotran/src/pflotran/pflotran
```

You should see an ~80–200 MB binary.

### 3.3 Verify baseline build

Run an existing PFLOTRAN regression test:

```bash
cd ~/pflotran/regression_tests/general
~/pflotran/src/pflotran/pflotran -input_prefix 1d_flux 2>&1 | tail -10
```

You should see at the end:

```
PMCSubsurfaceFlow
 Total Time: ... seconds
 ...
 Wall Clock Time: ... [sec]
```

If you see this, PFLOTRAN baseline works. Move on to Step 4.

### Troubleshooting Step 3 failures

**Symptom: "Cannot open module file 'hdf5.mod'"**

→ PETSc was built without `--download-hdf5-fortran-bindings=yes`.
   Re-do Step 2 with that flag.

**Symptom: `preconditioner_cpr.F90` errors about MPI_Comm or
c_int32_t**

→ PETSc version mismatch. Verify with `cd ~/petsc && git describe --tags`.
   It must show `v3.24.5`. If not, re-do Step 2.

**Symptom: build succeeds but pflotran binary is missing**

→ Linker error. Look at the last 50 lines of
   `/tmp/pflotran_baseline_build.log` for "undefined reference"
   lines. These usually indicate a missing PETSc subsystem.

---

## Step 4 — Get the EOS code

Time: 5 minutes.

### 4.1 Clone the EOS repository

```bash
cd ~
git clone https://github.com/Yongqiang100/pflotran-air-eos.git
cd pflotran-air-eos
ls
```

You should see:

```
CHANGELOG.md  CITATION.cff  docs  LICENSE
pflotran_integration  README.md  regression_tests  src
```

### 4.2 Verify the repository structure

```bash
find ~/pflotran-air-eos -type f -not -path '*/.git/*' | sort
```

You should see all the expected files (modules, tests, documentation).

---

## Step 5 — Verify standalone modules

Time: 10 minutes.

These tests confirm the EOS physics is correct *before* you
integrate it into PFLOTRAN. If a standalone test fails, the
problem is in the physics modules, not the integration.

### 5.1 Build the standalone tests

```bash
cd ~/pflotran-air-eos/src
make
ls -l test_*
```

You should see six executable files:

```
test_air_eos_module
test_air_eos_pr
test_air_eos_pr_sweep
test_air_henry
test_air_properties
test_air_watervapour
```

### 5.2 Run the tests

```bash
cd ~/pflotran-air-eos/src
for t in test_air_eos_pr test_air_henry test_air_watervapour \
         test_air_properties test_air_eos_module; do
  echo "=== $t ==="
  ./$t 2>&1 | tail -5
done
```

Each test should report `PASS` or similar success indicator at the
end. Take note of any failures — these are bugs to fix before
integration.

### 5.3 Run the cross-validation (optional)

The cross-validation compares the Fortran sweep output to an
independent Python implementation of PR + Peneloux.

**Important**: `test_air_eos_pr_sweep` writes `air_eos_sweep.csv`
itself, directly to the current directory. Do **not** redirect
stdout to that same filename — the redirect collides with the
program's own file write and corrupts the CSV. Just run the
program:

```bash
cd ~/pflotran-air-eos/src
./test_air_eos_pr_sweep
ls -lh air_eos_sweep.csv          # expect ~70 KB, 441 lines
wc -l air_eos_sweep.csv           # expect 441 (1 header + 440 data)
python3 cross_validate.py air_eos_sweep.csv
```

Expected output ends with:

```
Compared 440 grid points
  Max rel err rho   : 7.3e-08
  Max rel err Z     : 5.0e-08
  Max rel err phi   : 4.9e-08
  Max rel err h_dep : 6.1e-08
PASS: agreement to within 1e-6 (acceptable; differences are at
      floating-point roundoff scale)
```

The errors should be at 1e-7 level (essentially floating-point
roundoff between Fortran and Python). Any value above 1e-6
suggests a real numerical discrepancy and the cross_validate.py
script will report FAIL.

If the script reports `Compared 0 grid points` followed by
"ERROR", the CSV file is corrupted. The most common cause is a
stale `air_eos_sweep.csv` from a previous redirect-collision run.
Delete it and re-run:

```bash
rm -f air_eos_sweep.csv
./test_air_eos_pr_sweep
python3 cross_validate.py air_eos_sweep.csv
```

---

## Step 6 — Integrate EOS modules into PFLOTRAN

Time: 30 minutes.

### 6.1 Copy modules into PFLOTRAN source

```bash
cd ~/pflotran-air-eos/src
cp air_*.F90 ~/pflotran/src/pflotran/
ls ~/pflotran/src/pflotran/air_*.F90
```

You should see five files: `air_eos_module.F90`,
`air_eos_pr_module.F90`, `air_henry_module.F90`,
`air_properties_module.F90`, `air_watervapour_module.F90`.

### 6.2 Update PFLOTRAN's object file list

```bash
cd ~/pflotran/src/pflotran
cp pflotran_object_files.txt pflotran_object_files.txt.bak
```

Edit `pflotran_object_files.txt`:

```bash
nano pflotran_object_files.txt
```

Find the line that starts:

```
eos_obj = \
```

Look at the existing entries. They follow this pattern:

```
eos_obj = \
        ${common_src}co2_sw.o \
        ${common_src}co2_span_wagner_spline.o \
        ${common_src}eos.o \
        ${common_src}eos_gas.o \
        ${common_src}gas_eos_mod.o
```

Add the five air entries in alphabetical order, before the
existing entries. The complete block should look like:

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

Save (Ctrl+O, Enter, Ctrl+X).

Verify:

```bash
grep "air_" pflotran_object_files.txt
```

Should show five lines for `air_*.o`.

### 6.3 Update PFLOTRAN's dependency tracking

```bash
cp pflotran_dependencies.txt pflotran_dependencies.txt.bak
```

Edit `pflotran_dependencies.txt`:

```bash
nano pflotran_dependencies.txt
```

Add at the end:

```
air_eos_module.o : \
  air_eos_pr_module.o \
  air_henry_module.o \
  air_properties_module.o \
  air_watervapour_module.o
```

Save and verify:

```bash
grep -A 5 "air_eos_module" pflotran_dependencies.txt
```

### 6.4 Apply patches to eos_gas.F90 and eos.F90

The patches are documented in
`~/pflotran-air-eos/pflotran_integration/integration_patch.txt`.
Make a backup of the files before editing:

```bash
cd ~/pflotran/src/pflotran
cp eos_gas.F90 eos_gas.F90.bak
cp eos.F90 eos.F90.bak
```

Open the patch file in a viewer:

```bash
less ~/pflotran-air-eos/pflotran_integration/integration_patch.txt
```

Read it carefully. There are 5 sections:

- **A1**: Add 2 `use` statements near the top of `eos_gas.F90`
- **A2**: Add 1 entry to the public list in `eos_gas.F90`
- **A3**: Add the setter routine to `eos_gas.F90`
- **A4**: Add 4 adapter routines to `eos_gas.F90`
- **B1**: Add 2 lines to the keyword dispatcher in `eos.F90`

Apply each section by editing the relevant file with `nano` or
your preferred editor. Be careful about indentation: match the
surrounding code's indentation exactly.

### 6.5 Verify edits

After applying all patches:

```bash
cd ~/pflotran/src/pflotran

# Should find the new use statements (A1)
grep -n "use Air_EOS_module" eos_gas.F90
grep -n "use Air_Properties_module" eos_gas.F90

# Should find the new setter in the public list (A2)
grep -n "EOSGasSetPRPenelouxAir" eos_gas.F90

# Should find the new keyword case (B1)
grep -n "PR_PENELOUX_AIR" eos.F90
```

Each `grep` should return one or more matches. If any are empty,
that patch section wasn't applied correctly.

---

## Step 7 — Build PFLOTRAN with EOS modules

Time: 30–60 minutes.

```bash
cd ~/pflotran/src/pflotran
make clean
make pflotran 2>&1 | tee /tmp/pflotran_eos_build.log | tail -30
```

Successful indicators:

- The `air_*.o` object files compile successfully early in the build
- No errors during compilation of `eos_gas.F90` or `eos.F90`
- Final linking produces a `pflotran` executable

If you see compilation errors specific to the air modules, double-
check the patch application in Step 6.

Common errors and fixes:

**"Parameter 'dp' has not been declared"**

→ A local variable named `dP` (uppercase P) collided with a kind
   parameter `dp`. Rename local variables to `dP_step` or similar
   in the adapter routines.

**"Symbol 'X' has no IMPLICIT type"**

→ A `use` statement is missing. Verify all `use` statements in
   patch section A1 are applied.

**Linker error "undefined reference to air_..."**

→ The `air_*.o` files weren't added to `pflotran_object_files.txt`.
   Re-do Step 6.2.

Verify the new binary works:

```bash
~/pflotran/src/pflotran/pflotran -version 2>&1 | head -5
```

Should print version info without crashing.

Re-run the baseline test to confirm no regression:

```bash
cd ~/pflotran/regression_tests/general
~/pflotran/src/pflotran/pflotran -input_prefix 1d_flux 2>&1 | tail -10
```

Same wall clock time as in Step 3.3.

---

## Step 8 — Numerical validation

Time: 5 minutes.

Create the validation deck:

```bash
mkdir -p ~/pflotran/regression_tests/general/air_caes
cd ~/pflotran/regression_tests/general/air_caes
```

Copy the validation deck from the repository:

```bash
cp ~/pflotran-air-eos/regression_tests/validate_pr_air.in .
ls *.in
```

Run it:

```bash
~/pflotran/src/pflotran/pflotran -input_prefix validate_pr_air 2>&1 | tail -10
```

Should see "Wall Clock Time" at the end. Then check the output:

```bash
cat validate_pr_air-001.tec
```

You should see one data line. Check the values match the reference:

| Column | Field | Reference value | Tolerance |
|--------|-------|-----------------|-----------|
| 4 | Gas Density [kg/m³] | 87.92 | ±0.5 |
| 5 | Gas Viscosity [Pa·s] | 2.063e-5 | ±5e-7 |
| 6 | Gas Pressure [Pa] | 8.000e6 | exact |
| 7 | Gas Saturation | 0.95 | exact |
| 8 | Temperature [°C] | 40.0 | exact |

If gas density is ~87.92, the integration is numerically correct.

If gas density is off by a factor of ~30, you have a kg/m³ vs
kmol/m³ unit-conversion bug in the adapter routines. Review
Step 6.4 A4 carefully.

---

## Step 9 — 1D flow simulation

Time: 5 minutes.

This tests the EOS under real flow conditions, not just static
single-cell evaluation.

```bash
cd ~/pflotran/regression_tests/general/air_caes
cp ~/pflotran-air-eos/regression_tests/level3_1d_injection.in .
~/pflotran/src/pflotran/pflotran -input_prefix level3_1d_injection 2>&1 | tail -30
```

Look for:

- 70 time steps completed (`Step 70`)
- ~155 newton iterations total (about 2 per step — quadratic convergence)
- 0 timestep cuts (`cuts = 0`)
- Wall clock < 1 second

This indicates the analytic derivatives in the EOS are correct.

Examine the gas saturation profile:

```bash
head -10 level3_1d_injection-004.tec
```

The first column is position (x in meters). The gas saturation
column (7th value on each data line) should show:

- Near x = 0.5 m: saturation > 0.6 (gas front)
- Near x = 4.5 m: saturation ~ 0.05 (back to initial)
- Beyond x ~ 5 m: saturation ~ 0.05 (unperturbed)

This confirms the gas front advanced ~4 m in 1 hour, which is
physically correct for a 1000 mD aquifer with 5 bar driving
pressure.

---

## Completion checklist

Tick each item as you confirm:

- [ ] Step 1: gfortran 11+, gcc 11+, git, make, python3 installed
- [ ] Step 2.1: PETSc cloned and checked out at v3.24.5
  (verify: `cd ~/petsc && git describe --tags`)
- [ ] Step 2.3: PETSc built (`libpetsc.*` exists in arch-linux-c-opt/lib/)
- [ ] Step 2.4: `make check` shows "Completed test examples"
- [ ] Step 2.5: PETSC_DIR and PETSC_ARCH environment variables set
- [ ] Step 3.2: PFLOTRAN baseline binary exists
- [ ] Step 3.3: 1d_flux regression test passes
- [ ] Step 4: pflotran-air-eos repo cloned
- [ ] Step 5: All 6 standalone tests pass
- [ ] Step 6: All 5 air modules copied, all patches applied
- [ ] Step 7: PFLOTRAN with EOS rebuilds without errors
- [ ] Step 7: 1d_flux regression test still passes
- [ ] Step 8: validate_pr_air shows gas density ~87.92 kg/m³
- [ ] Step 9: level3_1d_injection shows ~70 steps, ~2 newton/step, 0 cuts

When all boxes are ticked, the reproduction is complete.

---

## Total disk usage at end of reproduction

```bash
du -sh ~/petsc ~/pflotran ~/pflotran-air-eos
df -h ~
```

Expected:

- `~/petsc`: ~3 GB
- `~/pflotran`: ~1 GB
- `~/pflotran-air-eos`: ~25 MB
- Total: ~4 GB

If you need to free space, the largest deletable item is
`~/petsc/arch-linux-c-opt/externalpackages/` (~2 GB of build
artifacts that aren't needed after make completes).

---

## Where to get help

If you get stuck at any step:

1. Re-read the troubleshooting notes for that step in this document
2. Check PFLOTRAN's own troubleshooting page:
   https://www.pflotran.org/documentation/user_guide/how_to/faq.html
3. Search the PFLOTRAN users mailing list:
   https://groups.google.com/g/pflotran-users
4. Open an issue at the pflotran-air-eos repository:
   https://github.com/Yongqiang100/pflotran-air-eos/issues

---

## Version of this document

This guide was written and tested against:

- PETSc v3.24.5 (released early 2026)
- PFLOTRAN main branch as of May 2026
- pflotran-air-eos v0.1.x

If PFLOTRAN's official documentation now recommends a different
PETSc version, use that instead — the structure of these steps
remains valid.
