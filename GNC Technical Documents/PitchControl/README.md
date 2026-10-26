# Pitch-control moment analysis

This MATLAB workflow prepares one trajectory/aerodynamic history and calculates vehicle, static-trim, and dynamic control moments for configurable dimensional CG locations and pitch inertias. Original input data, `InterpolateTraj.m`, and `interpolateAeroTrajectory.m` are preserved. The existing aerodynamic lookup helper is reused directly from the parent folder.

## Run

From the project root in MATLAB:

```matlab
addpath(fullfile('GNC Technical Documents','PitchControl'));
results = runPitchControlExample();
```

Edit **`pitchControlConfig.m`** for paths, CGs, inertias, normalization references, interpolation settings, and source confirmations. Each run creates a new `results/run_<timestamp>/` directory, so earlier studies and source files are not overwritten. MATLAB R2024a and Curve Fitting Toolbox (`csaps`, `fnval`, `fnder`) are the tested target. The exporter also uses MATLAB's bundled JVM for XML styling. No Excel installation or Python/JavaScript dependency is required to run the workflow.

## Inputs and reference status

The configuration defaults to the new Bulge body aero CSV and `bulge_30km_max_collocation.csv`. `Sref_m2` is read from the aero CSV (currently 0.02555495 m^2), as is mass (currently 120 kg). Mass is metadata only; it does not establish the optimizer's mass or determine `Iyy`.

The example uses `xCG_m = [1.08 1.12 1.16 1.20 1.24]`, `Iyy_kg_m2 = 15.536`, `Cbar_m = 0.18038`, and `xMomentReference_m = 0`. These are **provisional inherited values**, not independently confirmed properties of the new geometry. The legacy script instead uses Cbar = 0.18188 m. The new CSV supplies no Cbar or independently confirmed moment-reference/sign metadata.

Sref and Cbar must match the coefficient normalization, not just an arbitrary physical area or chord. An Sref override that disagrees with the database is rejected. Cbar remains an explicit input whose source requires confirmation. No automatic coefficient renormalization is performed. Set `config.verification.<field>.source` and `.confirmed` after actual engineering confirmation; numerical tests do not confirm physical conventions. Fields are `Sref`, `Cbar`, `Iyy`, `CG`, `MomentReference`, `CoordinateDirection`, and `PitchSign`.

The implemented moment translation assumes aft-positive x and zero vertical offset. Confirmation records document these assumptions; changing a note does not implement a different sign convention. Inertia is positive, constant in time, and about the pitch axis through the selected CG. No inertia is inferred from mass or shifted automatically as CG changes.

## Change cases without repeating interpolation

```matlab
config = pitchControlConfig();
prepared = preparePitchTrajectory(config);

% Pair each CG with an explicitly supplied inertia.
config.xCG_m = [1.08 1.16 1.24];
config.Iyy_kg_m2 = [14 15.536 17]; % Illustrative, not measured inertias.
config.caseMode = "paired";
config.caseLabels = ["Forward" "Middle" "Aft"];
paired = evaluatePitchCases(prepared, config);

% Evaluate all 3 x 3 combinations without preparing the trajectory again.
config.caseMode = "cartesian";
config.caseLabels = strings(0,1); % Generate unique IDs for all nine cases.
grid = evaluatePitchCases(prepared, config);

% Full exports for these cases; skip repeated sensitivity work if desired.
config.runSensitivity = false;
results = calculatePitchControlMoments(config, prepared);
```

A scalar Iyy is broadcast to every CG in paired mode. A paired vector must have the same length as xCG. Cartesian ordering varies xCG fastest, then Iyy. IDs never round coordinates. All cases share exactly the same common trajectory and aerodynamic columns.

`evaluatePitchCases` only evaluates cases in memory: it never interpolates or writes files. `calculatePitchControlMoments` adds optional sensitivity studies and exports. Its prepared-history reuse guard rejects changed filenames, source sizes/modification times, dt, or smoothing settings; regenerate preparation when those change. `evaluatePitchCases` deliberately accepts an in-memory prepared snapshot, including an archived one, without rereading source files. A snapshot cannot establish that an external file has remained unchanged.

## Interpolation

- Preserve the raw trajectory and all original columns; time-sort a separate copy and reject duplicate times.
- Build an endpoint-inclusive uniform grid with `ceil(duration/dt)` intervals and actual spacing no larger than dt (default 0.25 s).
- Convert alpha and gamma to radians, fit separate cubic smoothing splines (`smoothParam=0.999`), then evaluate their first and second analytical derivatives. The planar approximation is theta = alpha + gamma.
- PCHIP-interpolate Mach, altitude, speed, and trajectory dynamic pressure. Convert kPa to Pa exactly once.
- Query the existing linear aerodynamic helper at paired `(Mach, altitude, smoothed alpha)` states. Validate the complete grid and prohibit extrapolation; the helper only tolerates boundary roundoff of `1e-10*max(1, axis span)`.
- Keep database CL/CD as `CL`/`CD`; keep PCHIP versions of original trajectory coefficients as `TrajectoryCL_pchip`/`TrajectoryCD_pchip`. At original unsmoothed states, compare supplied trajectory CL/CD against database lookup in `Node_Diagnostics`.
- Calculate LD from interpolated CL/CD. Neither database `q_kPa` nor `xcp_over_L` is used to reconstruct the moment.

Sensitivity diagnostics compare the baseline, half timestep, and configured alternative smoothing parameters (default 0.99 and 1). Each variant interpolates once for all cases. It does not retune, clip peaks, or replace the baseline. Unavailable variants are explicitly recorded with reasons and blank numerical results. A relative peak change is also blank if the baseline peak is zero. Dense sampling does not recover dynamics absent from the original collocation nodes.

## Equations and scope

```text
Cm_CG = Cm_nose + CN*(xCG - xMomentReference)/Cbar
K = q_Pa*Sref*Cbar
M_vehicle = K*Cm_CG
M_inertial = Iyy*pitchAccel
M_trim_required = -M_vehicle
M_maneuver_required = M_inertial - M_vehicle
DeltaCm_trim = -Cm_CG
DeltaCm_maneuver = M_maneuver_required/K
```

The maneuver requirement already includes the trim contribution. Static trim is independent of inertia. Positive and negative mean the project's algebraic pitch directions, not an unconfirmed nose-up/down label. The model includes only the supplied clean aerodynamic moment and the prescribed planar acceleration. Vertical reference offsets, damping, thrust moments, six-degree-of-freedom coupling, time-varying inertia, and trajectory reoptimization are not modeled.

## Outputs

- `common_trajectory_aero.csv`: one common state/kinematic/aerodynamic history.
- `case_history.csv`: long-format histories for every case, including clean vehicle, inertial, static, and dynamic control moments and coefficients.
- `summary.csv`: independent minimum, maximum, and maximum-absolute values for each moment/coefficient, signed value, magnitude, sign-presence flags, time/history indices, and flight conditions. Ties select the earliest row. A minimum may be positive if no negative values occur; consult `HasNegative`/`HasPositive`.
- `node_diagnostics.csv`: original-state CL/CD comparison and alpha/gamma fit error.
- `interpolation_sensitivity.csv`: per-case peak changes under fixed diagnostic variants, when enabled.
- `quality_checks.csv`: numerical checks and distinct physical-reference confirmation records.
- `pitch_control_results.xlsx`: required input, raw, common, case, summary, and quality sheets, plus cases, node diagnostics, aerodynamic grid, sensitivity, and an editable formula example.
- `pitch_control_results.mat`: `results` structure, including reusable preparation, configuration, cases, histories, diagnostics, and output paths.
- `plots/`: PNG and MATLAB FIG versions of vehicle/static/dynamic moments, control coefficients versus Mach, pitch motion, spline fits, and sensitivity. Plot sample order remains chronological, including nonmonotonic Mach.
- Per-run `README.md`: settings, sources, reference status, and interpretation.

Main workbook sheets are **frozen numerical MATLAB exports**. Editing them does not update the calculation. Equations and units are recorded in `Inputs_Verification`. `Formula_Example` contains three sample states for the first case, with genuine Excel formulas in columns K:Q, cached MATLAB values, and automatic recalculation on open. Those examples are independent of the main summary. The runner verifies saved table dimensions and case-history values, plus formula insertion; it does not claim to execute native Excel recalculation.

Vehicle moment extrema are encountered clean aerodynamic moments, not available control authority. Required control moments and their coefficients must later be checked against a control-surface model at each flight condition, including both signs. This implementation does not calculate deflection, actuator loads, or prove trim feasibility.

## Tests

```matlab
addpath(fullfile('GNC Technical Documents','PitchControl'));
testResults = runPitchControlTests();
```

Tests cover dimensional CG and inertia sensitivities, zero acceleration, static inertia independence, paired/Cartesian behavior, raw/source units, a manufactured multi-affine aerodynamic field, source-node lookup with reordered rows, out-of-range and invalid grid errors, independent extrema, absent signs, reference mismatches, preparation reuse, duplicate/unsorted time, and small workbook/CSV exports including executable formula XML. Tests write only disposable fixtures and new temporary export directories. The original project files are not modified.
