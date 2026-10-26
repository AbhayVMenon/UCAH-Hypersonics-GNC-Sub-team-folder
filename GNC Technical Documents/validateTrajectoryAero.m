function validateTrajectoryAero()
% Validate grid indexing, paired queries, rejection behavior and CG shift.
folder = fileparts(mfilename('fullpath'));
aero = readtable(fullfile(folder, '..', 'Designs', 'Lifting Body V2', ...
    'Aero coeffs.csv'));
names = {'CN', 'CA', 'CL', 'CD', 'CDpres', 'CDfric', 'Cm_nose'};

% Recover every original node, with source rows deliberately reversed.
reversed = aero(end:-1:1, :);
atNodes = interpolateAeroTrajectory(reversed, aero.Mach, aero.alt_m, aero.alpha_deg);
assert(max(abs(atNodes{:, names} - aero{:, names}), [], 'all') < 1e-11);

% Independent trilinear midpoint check: all eight corners have equal weight.
corners = ismember(aero.Mach, [4, 5]) & ...
    ismember(aero.alt_m, [15000, 20000]) & ismember(aero.alpha_deg, [4, 5]);
assert(nnz(corners) == 8);
middle = interpolateAeroTrajectory(aero, 4.5, 17500, 4.5);
assert(max(abs(middle{:, names} - mean(aero{corners, names}, 1))) < 1e-11);

% Row inputs still represent two paired queries and return two output rows.
pair = interpolateAeroTrajectory(aero, [4, 5], [15000, 20000], [4, 5]);
assert(isequal(size(pair), [2, 8]));
mustFail(@() interpolateAeroTrajectory(aero, 2.9, 17500, 4.5), 'Aero:OutOfRange');
mustFail(@() interpolateAeroTrajectory(aero, 4.5, 37000, 4.5), 'Aero:OutOfRange');
mustFail(@() interpolateAeroTrajectory(aero, 4.5, 17500, 26), 'Aero:OutOfRange');
mustFail(@() interpolateAeroTrajectory(aero(2:end,:), 4.5, 17500, 4.5), ...
    'Aero:IncompleteGrid');
mustFail(@() interpolateAeroTrajectory([aero; aero(1,:)], 4.5, 17500, 4.5), ...
    'Aero:DuplicateState');

out = readtable(fullfile(folder, 'pitch_maneuver_history_with_aero.csv'));
assert(width(out) == 24 && all(isfinite(out{:,:}), 'all'));
assert(all(diff(out.Time_s) > 0));
cg = [0.50, 0.54, 0.58, 0.62, 0.66, 0.70];
for j = 1:numel(cg)
    % Interpolating shifted source moments must match shifting interpolated moments.
    shifted = aero;
    shifted.Cm_nose = aero.Cm_nose + aero.CN .* (2.0164 * cg(j) / 0.18188);
    shiftedResult = interpolateAeroTrajectory(shifted, out.Mach, ...
        out.Altitude_m, out.Alpha_deg);
    name = sprintf('CmClean_h%03d', round(100 * cg(j)));
    assert(max(abs(out.(name) - shiftedResult.Cm_nose)) < 1e-10);
end
fprintf(['PASS: all source nodes, shuffled input, 8-corner midpoint, paired queries, ' ...
    'out-of-range/missing/duplicate rejection, all six CG shifts.\n']);
fprintf('Output: %d rows x %d columns, all finite.\n', height(out), width(out));
end

function mustFail(f, identifier)
try
    f();
catch problem
    assert(strcmp(problem.identifier, identifier), ...
        'Expected %s; got %s.', identifier, problem.identifier);
    return;
end
error('Expected failure: %s', identifier);
end
