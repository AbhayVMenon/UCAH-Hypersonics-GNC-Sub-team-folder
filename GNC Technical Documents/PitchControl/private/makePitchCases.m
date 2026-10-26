function cases = makePitchCases(config)
%MAKEPITCHCASES Case IDs never encode rounded geometry or inertia.
if config.caseMode == "cartesian"
    [x, inertia] = ndgrid(config.xCG_m, config.Iyy_kg_m2);
    x = x(:);
    inertia = inertia(:);
elseif isscalar(config.Iyy_kg_m2)
    x = config.xCG_m;
    inertia = repmat(config.Iyy_kg_m2, numel(x), 1);
else
    x = config.xCG_m;
    inertia = config.Iyy_kg_m2;
end
n = numel(x);
ids = compose("Case_%04d", (1:n).');
labels = config.caseLabels;
if isempty(labels)
    labels = ids;
end
assert(numel(labels) == n && numel(unique(labels)) == n, 'Pitch:Labels', ...
    'Supply exactly %d unique case labels or leave caseLabels empty.', n);
cases = table(ids, labels, x, inertia, 'VariableNames', ...
    {'CaseID','CaseLabel','xCG_m','Iyy_kg_m2'});
end
