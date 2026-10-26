function config = pitchControlConfig()
%PITCHCONTROLCONFIG Editable starting configuration for the bulge-body study.
% Values inherited from the handoff are PROVISIONAL, not new-geometry data.
projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
config.trajectoryFile = fullfile(projectRoot, 'Designs', 'Bulge body', ...
    'Trajectory data', 'bulge_30km_max_collocation.csv');
config.aeroFile = fullfile(projectRoot, 'Designs', 'Bulge body', ...
    'MATLAB Data', 'cst_v7_bulge_nxV1p1_otis_deck.csv');
config.outputDirectory = fullfile(fileparts(mfilename('fullpath')), 'results');

%% Design cases: dimensional xCG in metres; Iyy about each case's CG.
config.xCG_m = [1.08 1.12 1.16 1.20 1.24];
config.Iyy_kg_m2 = 15.536;
config.caseMode = "paired";
config.caseLabels = strings(0, 1);

%% Coefficient reference geometry: must match the original aero normalization.
config.Sref_m2 = [];          % Read from the NEW aero CSV, never the old deck.
config.Cbar_m = 0.18038;      % PROVISIONAL: old script instead uses 0.18188 m.
config.xMomentReference_m = 0; % PROVISIONAL: assumed sharp-nose reference.

%% Interpolation and diagnostics.
config.dt_s = 0.25;
config.smoothParam = 0.999;
config.runSensitivity = true;
config.sensitivitySmoothParams = [0.99 1.0]; % Diagnostic, not an optimizer.
config.coefficientComparisonTolerance = 1e-4; % Absolute CL/CD diagnostic only.
config.exportResults = true;
config.createPlots = true;
config.maxCasesPerPlot = 6;

%% Physical confirmation is manual and independent of numerical checks.
names = {'Sref','Cbar','Iyy','CG','MomentReference','CoordinateDirection','PitchSign'};
sources = { ...
    'Sref_m2 in supplied aero CSV; dataset-owner confirmation pending', ...
    'Handoff: 0.18038 m; old InterpolateTraj.m: 0.18188 m; deck confirmation pending', ...
    'Handoff: 15.536 kg m^2; confirm pitch inertia about each selected CG', ...
    'Handoff candidates in metres; confirm origin and new-geometry applicability', ...
    'Cm_nose column name; assumed xref=0 is not independent confirmation', ...
    'Provisional aft-positive x; longitudinal shift only, zero vertical offset', ...
    'Positive/negative pitch only; physical nose-up/nose-down mapping unconfirmed'};
for k = 1:numel(names)
    config.verification.(names{k}) = struct('source', sources{k}, 'confirmed', false);
end
end
