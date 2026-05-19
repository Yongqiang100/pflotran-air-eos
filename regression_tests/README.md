# `regression_tests/` — PFLOTRAN input decks

Test decks demonstrating the high-pressure air EOS at increasing
physical complexity. Run these after completing the PFLOTRAN
integration (see [`../docs/WORKFLOW.md`](../docs/WORKFLOW.md)).

## Test ladder

| Level | Geometry | Physics | Expected runtime |
|---|---|---|---|
| 1 | Single cell | Static, isothermal | seconds |
| 2 | Single cell | Pressure transient, isothermal | seconds |
| 3 | 1D column, 100 m | Single-pass gas injection | < 1 minute |
| 5 | 2D axisymmetric (Yang geometry) | One daily cycle, isothermal | < 10 minutes |

## How to run

After PFLOTRAN integration is complete:

```bash
$PFLOTRAN_DIR/src/pflotran/pflotran -input_prefix level1_single_cell
$PFLOTRAN_DIR/src/pflotran/pflotran -input_prefix level2_pressure_relax
$PFLOTRAN_DIR/src/pflotran/pflotran -input_prefix level3_1d_injection
$PFLOTRAN_DIR/src/pflotran/pflotran -input_prefix level5_2d_yang_singlecycle
```

Each level produces output files (`.tec`, `.h5`, or `.out` depending
on the OUTPUT block format). Expected results documented at the
bottom of each `.in` file.

## Verification values

Standalone wrapper output at 80 bar / 40 °C / 0.5 mol/kg NaCl:

```
Gas density       = 87.95 kg/m³
Gas viscosity     = 2.063 × 10⁻⁵ Pa·s
Fugacity coeff    = 1.003
Water vapour y_w  = 9.55 × 10⁻⁴
Dissolved O₂      = 14.57 mol/m³ (466 mg/L)
```

PFLOTRAN should reproduce these values to better than 0.1 % when the
deck is set up correctly.

## Levels not in this initial release

The full Yang et al. (2024) reproducibility (level 6 = 100 cycles)
and geochemistry-coupled simulation (level 7 = full reactive
transport) are not included in this initial release. Outlines for
both are in the project READMEs; flesh them out once levels 1-5
are confirmed working.

## Important notes

- These decks use the `PR_PENELOUX_AIR` keyword, which requires the
  integration patches in [`../pflotran_integration/`](../pflotran_integration/)
  to be applied first
- Some decks were adapted from the initial deliverable-6 versions to
  match the actual integration approach (function-pointer dispatch via
  `EOS GAS / DENSITY PR_PENELOUX_AIR` rather than a new fluid type)
- The 2D axisymmetric level 5 deck uses a simplified flat geometry
  rather than Yang et al.'s tilted parallelogram. For full Yang
  reproduction, add the dip via rotated gravity vector
