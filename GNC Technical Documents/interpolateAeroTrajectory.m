function result = interpolateAeroTrajectory(aero, mach, altitude_m, alpha_deg)
%INTERPOLATEAEROTRAJECTORY Linear aero lookup at paired trajectory states.
% Inputs may be row or column vectors, but must have equal lengths.
% Requires a complete rectangular source grid. Extrapolation is prohibited.
% CL/CD are interpolated tabulated values; LD is calculated from their ratio.

axisNames = {'Mach', 'alt_m', 'alpha_deg'};
coefficientNames = {'CN', 'CA', 'CL', 'CD', 'CDpres', 'CDfric', 'Cm_nose'};
required = [axisNames, coefficientNames];
assert(all(ismember(required, aero.Properties.VariableNames)), ...
    'Aero:MissingColumns', 'The aero CSV is missing required columns.');
source = aero{:, required};
assert(isnumeric(source) && all(isfinite(source(:))), ...
    'Aero:InvalidSource', 'Required aero values must be finite numbers.');

queryInputs = {mach, altitude_m, alpha_deg};
for k = 1:3
    validateattributes(queryInputs{k}, {'numeric'}, ...
        {'vector', 'real', 'finite', 'nonempty'});
end
assert(numel(mach) == numel(altitude_m) && numel(mach) == numel(alpha_deg), ...
    'Aero:QuerySize', 'Mach, altitude and alpha must have equal lengths.');
query = [mach(:), altitude_m(:), alpha_deg(:)];

[machGrid, ~, iM] = unique(aero.Mach);
[altGrid, ~, iH] = unique(aero.alt_m);
[alphaGrid, ~, iA] = unique(aero.alpha_deg);
axes = {machGrid, altGrid, alphaGrid};
gridSize = cellfun(@numel, axes);
assert(all(gridSize >= 2), 'Aero:GridSize', ...
    'Each aero grid dimension must contain at least two values.');
sourceIndex = sub2ind(gridSize, iM, iH, iA);
assert(numel(unique(sourceIndex)) == height(aero), ...
    'Aero:DuplicateState', 'Duplicate (Mach, altitude, alpha) rows in aero CSV.');
assert(height(aero) == prod(gridSize), ...
    'Aero:IncompleteGrid', 'The aero CSV has missing grid combinations.');

for k = 1:3
    lo = axes{k}(1);
    hi = axes{k}(end);
    % Accept only floating-point roundoff at boundaries, not extrapolation.
    tol = 1e-10 * max(1, hi - lo);
    outside = find(query(:, k) < lo - tol | query(:, k) > hi + tol, 1);
    if ~isempty(outside)
        error('Aero:OutOfRange', ...
            'Trajectory row %d: %s = %.12g is outside [%.12g, %.12g].', ...
            outside, axisNames{k}, query(outside, k), lo, hi);
    end
    query(:, k) = min(hi, max(lo, query(:, k)));
end

values = zeros(size(query, 1), numel(coefficientNames));
for k = 1:numel(coefficientNames)
    gridValues = zeros(gridSize);
    gridValues(sourceIndex) = aero.(coefficientNames{k});
    F = griddedInterpolant(axes, gridValues, 'linear', 'none');
    % Column query vectors are paired points, not a Cartesian query grid.
    values(:, k) = reshape(F(query(:, 1), query(:, 2), query(:, 3)), [], 1);
end
assert(all(isfinite(values(:))), 'Aero:InvalidOutput', ...
    'The aerodynamic interpolation returned nonfinite values.');
result = array2table(values, 'VariableNames', coefficientNames);
assert(all(result.CD > 0), 'Aero:Drag', ...
    'Positive CD is required to calculate lift-to-drag ratio.');
result.LD = result.CL ./ result.CD;
end
