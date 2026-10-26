function config = validatePitchConfig(config)
%VALIDATEPITCHCONFIG Normalize public options without inventing reference data.
assert(isstruct(config) && isscalar(config), 'Pitch:Config', ...
    'config must be a scalar structure. Start with pitchControlConfig().');
required = {'xCG_m','Iyy_kg_m2','Cbar_m','xMomentReference_m'};
for k = 1:numel(required)
    assert(isfield(config, required{k}), 'Pitch:MissingConfig', ...
        'Missing required config field: %s.', required{k});
end
defaults = struct('caseMode', "paired", 'caseLabels', strings(0,1), ...
    'Sref_m2', [], 'dt_s', 0.25, 'smoothParam', 0.999, ...
    'runSensitivity', true, 'sensitivitySmoothParams', [0.99 1], ...
    'coefficientComparisonTolerance', 1e-4, 'exportResults', true, ...
    'createPlots', true, 'maxCasesPerPlot', 6, 'verification', struct());
names = fieldnames(defaults);
for k = 1:numel(names)
    if ~isfield(config, names{k})
        config.(names{k}) = defaults.(names{k});
    end
end
validateattributes(config.xCG_m, {'numeric'}, {'vector','real','finite','nonempty'});
validateattributes(config.Iyy_kg_m2, {'numeric'}, ...
    {'vector','real','finite','nonempty','positive'});
validateattributes(config.Cbar_m, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(config.xMomentReference_m, {'numeric'}, {'scalar','real','finite'});
validateattributes(config.dt_s, {'numeric'}, {'scalar','real','finite','>',0,'<=',0.25});
validateattributes(config.smoothParam, {'numeric'}, {'scalar','real','finite','>=',0,'<=',1});
validateattributes(config.maxCasesPerPlot, {'numeric'}, {'scalar','integer','positive','finite'});
validateattributes(config.coefficientComparisonTolerance, {'numeric'}, ...
    {'scalar','real','finite','nonnegative'});
validateattributes(config.sensitivitySmoothParams, {'numeric'}, ...
    {'vector','real','finite','nonempty','>=',0,'<=',1});
if ~isempty(config.Sref_m2)
    validateattributes(config.Sref_m2, {'numeric'}, {'scalar','real','finite','positive'});
end
for name = {'runSensitivity','exportResults','createPlots'}
    validateattributes(config.(name{1}), {'logical'}, {'scalar'});
end
config.caseMode = string(config.caseMode);
assert(isscalar(config.caseMode) && any(config.caseMode == ["paired","cartesian"]), ...
    'Pitch:CaseMode', 'caseMode must be "paired" or "cartesian".');
config.xCG_m = double(config.xCG_m(:));
config.Iyy_kg_m2 = double(config.Iyy_kg_m2(:));
if config.caseMode == "paired"
    assert(isscalar(config.Iyy_kg_m2) || numel(config.Iyy_kg_m2) == numel(config.xCG_m), ...
        'Pitch:PairedDimensions', 'Paired mode needs scalar Iyy or one Iyy per xCG.');
end
config.caseLabels = string(config.caseLabels(:));
assert(all(~ismissing(config.caseLabels) & strlength(strtrim(config.caseLabels)) > 0), ...
    'Pitch:Labels', 'Supplied case labels must be nonblank.');
assert(isstruct(config.verification) && isscalar(config.verification), ...
    'Pitch:Verification', 'verification must be a scalar structure.');
keys = {'Sref','Cbar','Iyy','CG','MomentReference','CoordinateDirection','PitchSign'};
for k = 1:numel(keys)
    key = keys{k};
    if ~isfield(config.verification, key)
        config.verification.(key) = struct('source', 'Not provided', 'confirmed', false);
    end
    entry = config.verification.(key);
    assert(isstruct(entry) && isscalar(entry) && ...
        isfield(entry,'source') && isfield(entry,'confirmed'), ...
        'Pitch:Verification', 'verification.%s needs source and confirmed.', key);
    validateattributes(entry.confirmed, {'logical'}, {'scalar'});
    entry.source = string(entry.source);
    assert(isscalar(entry.source) && ~ismissing(entry.source) && ...
        strlength(strtrim(entry.source)) > 0, 'Pitch:Verification', ...
        'verification.%s.source must be nonblank.', key);
    config.verification.(key) = entry;
end
end
