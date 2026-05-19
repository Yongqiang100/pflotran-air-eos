# `pflotran_integration/` — PFLOTRAN source patches

This directory contains the source-code patches needed to integrate
the standalone modules from [`../src/`](../src/) into PFLOTRAN's
GENERAL mode.

## Files

- `integration_patch.txt` — Documented patches with code blocks for:
  - 2 new `use` statements in `eos_gas.F90`
  - 1 entry in the public list of `eos_gas.F90`
  - 1 new setter routine `EOSGasSetPRPenelouxAir()`
  - 4 new adapter routines (`EOSGasDensityPRAir`, `EOSGasEnergyPRAir`,
    `EOSGasDensityEnergyPRAir`, `EOSGasViscosityPRAir`)
  - 1 new keyword case `'PR_PENELOUX_AIR'` in `eos.F90`

## How to apply

Full step-by-step procedure in [`../docs/WORKFLOW.md`](../docs/WORKFLOW.md)
Stages 3–8.

Briefly:

1. Build PFLOTRAN baseline and verify it works
2. Copy the five `.F90` files from `../src/` into `$PFLOTRAN_SRC/`
3. Add the new objects to `pflotran_object_files.txt` under `eos_obj`
4. Add module dependencies to `pflotran_dependencies.txt`
5. Apply the patches in `integration_patch.txt` to `eos_gas.F90` and
   `eos.F90`
6. Rebuild PFLOTRAN

## Activating in a deck

After successful integration:

```
EOS GAS
  DENSITY PR_PENELOUX_AIR
END
```

## Unit conventions in the adapters

PFLOTRAN's GENERAL mode internal units differ from the standalone
modules:

| Quantity | Standalone (SI) | PFLOTRAN GENERAL |
|---|---|---|
| Temperature | Kelvin | Celsius |
| Density | kg/m³ | kmol/m³ |
| Enthalpy | J/mol | J/kmol |
| Pressure | Pa | Pa |
| Viscosity | Pa·s | Pa·s |

The adapter routines in `integration_patch.txt` handle these
conversions explicitly. If you see PFLOTRAN reporting densities off
by ~30× or ~1000× from the standalone reference, suspect a
conversion bug in the adapters.

## Compatibility

These patches were developed against the PFLOTRAN main branch as of
May 2026. If you're using a different PFLOTRAN version, the line
numbers and surrounding context will differ; the *logical* patches
should still apply, but you'll need to read the live source to find
the right insertion points.

If the integration patches stop working with a future PFLOTRAN
version, the most likely culprits are:

1. New `auxvar` field names (the wrapper writes to PFLOTRAN's auxvar
   through the dummy interfaces, so PFLOTRAN-side changes there can
   affect us)
2. Function-pointer signature changes (the Dummy interfaces in
   `eos_gas.F90`)
3. Build-system reorganization (object lists or dependency tracking)

Updating the patches for a new PFLOTRAN version typically takes
1–2 hours of source reading once you've done it once.
