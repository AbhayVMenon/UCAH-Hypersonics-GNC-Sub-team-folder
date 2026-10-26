# Validation record

Validated on 2026-10-24 using MATLAB R2024a Update 6 with Curve Fitting Toolbox. The original data and legacy scripts were not modified.

## Executed checks

- `runPitchControlTests`: **18 passed, 0 failed** on the final regression run.
- MATLAB Code Analyzer: **0 messages across 15 new MATLAB files**.
- Actual-data example: 81 raw collocation points, 1,482 aerodynamic grid rows, 3,067 common time samples, 5 cases, 15,335 case-history rows, and 90 summary rows.
- Actual common time spacing: 0.249972158189 s; original endpoints retained.
- Original-state CL/CD consistency: maximum absolute differences 5.33e-15 and 7.22e-16, respectively.
- Spline fit errors at original nodes: maximum alpha error 0.0058226 deg and gamma error 0.000536164 deg.
- MATLAB CSV round trips checked shapes, IDs, and numerical values. Workbook checks covered all sheet dimensions and every case-history value.
- An independent Python evaluation checked all six moment/coefficient equations for all 15,335 exported case rows and all 90 summary-row references. Maximum scaled discrepancy was 6.84e-12 against a 1e-10 tolerance.
- A separate spreadsheet engine (Artifact Tool) imported the final workbook, verified 21 formula-example values, and confirmed the expected dynamic-moment change after changing sample Iyy. The workbook was not rewritten by that check. Native Microsoft Excel recalculation was not executed.
- Every workbook sheet was visually previewed. Plot files were inspected for labels, legends, and diagnostic content.

Example output directory: `results/run_20261024_012901_110/`.

## Provisional default results

These values use Cbar = 0.18038 m, Iyy = 15.536 kg m^2, xref = 0, and the new database Sref = 0.02555495 m^2. Values are rounded here; exports retain precision.

| xCG (m) | Signed dynamic moment at largest magnitude (N m) | Time (s) | Mach |
|---:|---:|---:|---:|
| 1.08 | -462.564 | 758.666 | 3.60414 |
| 1.12 | 343.270 | 766.415 | 3 |
| 1.16 | 423.176 | 766.415 | 3 |
| 1.20 | 503.081 | 766.415 | 3 |
| 1.24 | 582.987 | 766.415 | 3 |

These are encountered requirements, not available control authority. Both signs and their separate limiting conditions are retained in the full summary.

## Sensitivity and remaining engineering review

Halving the time step did not change the reported maximum-absolute dynamic moments for these five cases at exported precision. Changing smoothing from 0.999 to 0.99 or 1 changed those peak magnitudes by at most approximately 0.2793%. Peak pitch acceleration ranged from approximately 0.031515 to 0.037027 rad/s^2 across the tested smoothing settings, compared with 0.036179 rad/s^2 at baseline. The pitch-motion plot shows oscillations around sharp early and late trajectory changes; the sensitivity check does not certify those dynamics as physical.

The new CSV does not resolve the legacy Cbar discrepancy (0.18038 versus 0.18188 m) or confirm the moment origin/sign convention. CG coordinates and inertia must be checked for the new geometry, and inertia must be about each CG. Reference values remain explicitly provisional. No deflection, actuator sizing, available authority, or trim-feasibility conclusion has been calculated.
