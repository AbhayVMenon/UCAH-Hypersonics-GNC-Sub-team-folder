function results = calculatePitchControlMoments(config, prepared)
%CALCULATEPITCHCONTROLMOMENTS Interpolation, case analysis, checks, and exports.
% results = calculatePitchControlMoments(config)
% results = calculatePitchControlMoments(config, previousResults.prepared)
% For instant case-only exploration, use evaluatePitchCases(prepared,config).
config = validatePitchConfig(config);
if nargin < 2
    prepared = preparePitchTrajectory(config);
else
    validateReuse(prepared,config);
end
results = evaluatePitchCases(prepared,config);
fprintf('%s\n',results.status);
fprintf('Reference Cbar %.8g m (legacy conflict: 0.18038 vs 0.18188 m).\n',config.Cbar_m);
fprintf('Common history: %d states; %d cases; Sref %.10g m^2.\n', ...
    height(results.common),height(results.cases),results.metadata.Sref_m2);
if config.runSensitivity
    results.sensitivity = pitchControlSensitivity(prepared,config);
    failed = sum(results.sensitivity.Status ~= "OK");
    results.qualityChecks = [results.qualityChecks; pitchCheck('Interpolation sensitivity','REVIEW', ...
        sprintf('Inspect %d per-case diagnostic rows; %d unavailable variants. Baseline unchanged.',height(results.sensitivity),failed))];
else
    results.qualityChecks = [results.qualityChecks; pitchCheck('Interpolation sensitivity','NOT_RUN', ...
        'Disabled in config; preparation and case calculations still complete.')];
end
if config.exportResults
    results = exportPitchControlResults(results);
    fprintf('Results written to %s\n',results.outputDirectory);
end
end

function validateReuse(prepared,config)
assert(isstruct(prepared) && isfield(prepared,'preparationConfig') && ...
    isfield(prepared,'schemaVersion') && prepared.schemaVersion == 1, ...
    'Pitch:PreparedData','Unsupported prepared-history structure.');
old = prepared.preparationConfig;
for name = {'trajectoryFile','aeroFile','dt_s','smoothParam'}
    key = name{1};
    assert(isfield(config,key), 'Pitch:MissingConfig', 'Missing config.%s.',key);
    if endsWith(key,'File')
        same = strcmp(string(old.(key)),string(config.(key)));
    else
        same = isequal(old.(key),config.(key));
    end
    assert(same,'Pitch:StalePreparation', ...
        '%s changed. Re-run preparePitchTrajectory before reusing the history.',key);
end
fields = {'trajectory','aero'};
files = {config.trajectoryFile,config.aeroFile};
for k = 1:2
    info = dir(files{k});
    oldInfo = prepared.sourceInfo.(fields{k});
    assert(isscalar(info) && info.bytes == oldInfo.bytes && info.datenum == oldInfo.datenum, ...
        'Pitch:StalePreparation', 'Source file changed or is missing: %s.',files{k});
end
end
