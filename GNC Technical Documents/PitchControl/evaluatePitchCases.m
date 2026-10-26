function results = evaluatePitchCases(prepared, config)
%EVALUATEPITCHCASES Evaluate design cases without trajectory/aero interpolation.
config = validatePitchConfig(config);
assert(isstruct(prepared) && isfield(prepared,'common') && ...
    isfield(prepared,'metadata'), 'Pitch:PreparedData', ...
    'Use the structure returned by preparePitchTrajectory.');
common = prepared.common;
required = {'Time_s','Mach','Altitude_m','Alpha_deg','DynamicPressure_Pa', ...
    'PitchAccel_rad_s2','CN','Cm_nose'};
assert(istable(common) && all(ismember(required,common.Properties.VariableNames)), ...
    'Pitch:PreparedData', 'Prepared history is missing required columns.');
values = common{:,required};
assert(isnumeric(values) && isreal(values) && all(isfinite(values),'all') && ...
    height(common) >= 2 && all(diff(common.Time_s) > 0), ...
    'Pitch:PreparedData', 'Prepared history must be finite, real, and strictly increasing in time.');
assert(all(common.DynamicPressure_Pa > 0), 'Pitch:DynamicPressure', 'Dynamic pressure must be positive.');
sref = prepared.metadata.Sref_m2;
validateattributes(sref,{'numeric'},{'scalar','real','finite','positive'});
if ~isempty(config.Sref_m2)
    assert(abs(config.Sref_m2-sref) <= 1e-10*max(abs(sref),eps), ...
        'Pitch:SrefMismatch', 'Sref override disagrees with prepared aerodynamic normalization.');
end
config.Sref_m2 = sref;
cases = makePitchCases(config);
n = height(common);
K = common.DynamicPressure_Pa*sref*config.Cbar_m;
blocks = cell(height(cases),1);
for k = 1:height(cases)
    % Longitudinal x shift only, zero vertical offset, provisional aft-positive x.
    cm = common.Cm_nose + common.CN*((cases.xCG_m(k)-config.xMomentReference_m)/config.Cbar_m);
    vehicle = K.*cm;
    inertial = cases.Iyy_kg_m2(k)*common.PitchAccel_rad_s2;
    trim = -vehicle;
    maneuver = inertial-vehicle; % Includes trim; do not add trim a second time.
    calc = table(cm,vehicle,inertial,trim,maneuver,-cm,maneuver./K, ...
        'VariableNames',{'Cm_CG','M_vehicle_Nm','M_inertial_Nm', ...
        'M_trim_required_Nm','M_maneuver_required_Nm','DeltaCm_trim','DeltaCm_maneuver'});
    assert(all(isfinite(calc{:,:}),'all'), 'Pitch:NonfiniteResult', 'Nonfinite result in %s.',cases.CaseID(k));
    residual = vehicle+maneuver-inertial;
    assert(all(abs(residual) <= 1e-10*max(1,abs(vehicle)+abs(maneuver)+abs(inertial))), ...
        'Pitch:MomentBalance','Moment balance failed in %s.',cases.CaseID(k));
    blocks{k} = [repmat(cases(k,:),n,1),common,calc];
end
history = vertcat(blocks{:});
assert(height(history) == n*height(cases), 'Pitch:CaseRows', 'Case-history row count mismatch.');
results.config = config;
results.cases = cases;
results.common = common;
results.caseHistory = history;
results.summary = summarizePitchCases(history);
results.metadata = prepared.metadata;
results.rawTrajectory = prepared.rawTrajectory;
results.nodeDiagnostics = prepared.nodeDiagnostics;
results.aeroGrid = prepared.aeroGrid;
results.prepared = prepared;
results.qualityChecks = [prepared.qualityChecks; ...
    pitchCheck('Shared case history','PASS',sprintf('%d cases each use identical %d-row common history.',height(cases),n)); ...
    pitchCheck('Moment balance and finite results','PASS','M_vehicle + M_maneuver = Iyy*pitchAccel; all calculated moments finite.'); ...
    pitchCheck('Moment and coefficient extrema','PASS','Each quantity searched independently; signed values and original history indices retained.')];

keys = {'Sref','Cbar','Iyy','CG','MomentReference','CoordinateDirection','PitchSign'};
confirmed = true;
for k = 1:numel(keys)
    entry = config.verification.(keys{k});
    status = "UNCONFIRMED";
    if entry.confirmed, status = "USER_CONFIRMED"; end
    confirmed = confirmed && entry.confirmed;
    results.qualityChecks = [results.qualityChecks; ...
        pitchCheck("Physical reference: "+keys{k},status,entry.source)];
end
if confirmed
    results.status = "REFERENCES USER-CONFIRMED; REVIEW INTERPOLATION DIAGNOSTICS";
else
    results.status = "PROVISIONAL - UNVERIFIED REFERENCES - NOT FOR DESIGN RELEASE";
end
results.qualityChecks = [results.qualityChecks; pitchCheck('Legacy Cbar discrepancy','REVIEW', ...
    sprintf('Handoff 0.18038 m; legacy script 0.18188 m. Selected %.12g m; new CSV does not establish Cbar.',config.Cbar_m))];
results.sensitivity = table();
results.outputDirectory = "";
end
