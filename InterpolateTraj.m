clear; 
clc; 

%% Read trajectory
inputFile = fullfile( ...
"C:\Users\abhay\OneDrive - purdue.edu\(A) UCAH Hypersonics Comp\Designs\Lifting Body V2\bulge_results_collocation.csv"); %INPUT YOUR FILE PATH FOR OPTIMIZED TRAJECTORY

T = readtable(inputFile);
requiredColumns = {'Time_s', 'Mach', 'Altitude_m', 'Speed_mps', ...
    'DynamicPressure_kPa', 'FlightPath_deg', 'Alpha_deg'};
assert(all(ismember(requiredColumns, T.Properties.VariableNames)), ...
    'The trajectory CSV is missing required columns.');
requiredData = T{:, requiredColumns};
assert(isnumeric(requiredData) && all(isfinite(requiredData(:))), ...
    'Required trajectory columns must contain finite numeric values.');
assert(height(T) >= 3, 'At least three trajectory points are required.');
T = sortrows(T, "Time_s");

t = T.Time_s(:);

% Confirm that time is strictly increasing
if any(diff(t) <= 0)
    error("Trajectory time must be strictly increasing.");
end

%% Convert angles to radians
alphaRaw = deg2rad(T.Alpha_deg(:));
gammaRaw = deg2rad(T.FlightPath_deg(:));

%% Analysis settings
dt = 0.25;             % Initial uniform timestep, seconds
smoothParam = 0.999;   % Retained setting; physical smoothing needs validation

validateattributes(dt, {'numeric'}, {'scalar', 'real', 'finite', 'positive'});
validateattributes(smoothParam, {'numeric'}, ...
    {'scalar', 'real', 'finite', '>=', 0, '<=', 1});

% Include both endpoints while keeping a uniform spacing no greater than dt.
nIntervals = ceil((t(end) - t(1)) / dt);
tQuery = linspace(t(1), t(end), nIntervals + 1).';

%% Smooth alpha and gamma separately
alphaSpline = csaps(t, alphaRaw, smoothParam);
gammaSpline = csaps(t, gammaRaw, smoothParam);

%% Evaluate angles on uniform time grid
% reshape(..., [], 1) guarantees a column regardless of fnval orientation.
alpha = reshape(fnval(alphaSpline, tQuery), [], 1);
gamma = reshape(fnval(gammaSpline, tQuery), [], 1);

theta = alpha + gamma;

%% Analytical spline derivatives
alphaRate = reshape(fnval(fnder(alphaSpline, 1), tQuery), [], 1);
gammaRate = reshape(fnval(fnder(gammaSpline, 1), tQuery), [], 1);

alphaAccel = reshape(fnval(fnder(alphaSpline, 2), tQuery), [], 1);
gammaAccel = reshape(fnval(fnder(gammaSpline, 2), tQuery), [], 1);

% Body pitch rate and pitch acceleration
pitchRate = alphaRate + gammaRate;        % rad/s
pitchAccel = alphaAccel + gammaAccel;     % rad/s^2

%% Interpolate remaining trajectory variables
mach = interp1(t, T.Mach, tQuery, "pchip");
velocity = interp1(t, T.Speed_mps, tQuery, "pchip");

dynamicPressure = 1000 * interp1( ...
    t, T.DynamicPressure_kPa, tQuery, "pchip");  % Pa

altitude = interp1(t, T.Altitude_m, tQuery, "pchip");
mach = mach(:);
velocity = velocity(:);
dynamicPressure = dynamicPressure(:);
altitude = altitude(:);

%% Aerodynamic coefficients only at the trajectory flight conditions
aeroFile = fullfile(fileparts(inputFile), 'Aero coeffs.csv'); %INPUT YOUR FILE PATH FOR THE CBAERO COEFFICIENTS
aero = readtable(aeroFile);
% Each query is one (Mach, altitude in m, alpha in degrees) triplet.
% This does not generate a new sweep over every database altitude.
aeroTrajectory = interpolateAeroTrajectory(aero, mach, altitude, rad2deg(alpha));

% Preliminary reference settings from the current design discussion.
% Confirm CREF and the moment reference against the original CBAero deck.
L = 2.0164;                   % m, design-frame vehicle length, (INPUT FOR REQUIRED L, cref)
cref = 0.18188;               % m, inferred value, not yet deck-confirmed
xMomentReference_m = 0;       % ASSUMES Cm_nose is about the sharp nose
xcgCandidates = [0.50, 0.54, 0.58, 0.62, 0.66, 0.70]; % design-frame xcg/L
validateattributes(cref, {'numeric'}, {'scalar', 'finite', 'positive'});
validateattributes(L, {'numeric'}, {'scalar', 'finite', 'positive'});
xcg_m = L .* xcgCandidates;

% Longitudinal x-shift only; assumes no z-offset from moment reference.
% A nonzero z-offset also requires the axial-force moment contribution.
CmCleanTrajectory = aeroTrajectory.Cm_nose + aeroTrajectory.CN .* ...
    ((xcg_m - xMomentReference_m) ./ cref);
fprintf('Interpolated aero at %d trajectory states for %d CG candidates.\n', ...
    numel(tQuery), numel(xcgCandidates));

%% Check how far smoothing moved the original trajectory nodes
alphaAtNodes = reshape(fnval(alphaSpline, t), [], 1);
gammaAtNodes = reshape(fnval(gammaSpline, t), [], 1);

maxAlphaErrorDeg = max(abs(rad2deg(alphaAtNodes - alphaRaw)));
maxGammaErrorDeg = max(abs(rad2deg(gammaAtNodes - gammaRaw)));

fprintf("Maximum alpha fitting error: %.5f deg\n", maxAlphaErrorDeg);
fprintf("Maximum gamma fitting error: %.5f deg\n", maxGammaErrorDeg);

%% Validate output shapes and values before creating the table
outputColumns = {tQuery, mach, altitude, velocity, dynamicPressure, ...
    alpha, gamma, theta, pitchRate, pitchAccel};
N = numel(tQuery);
assert(all(cellfun(@(v) isequal(size(v), [N, 1]), outputColumns)), ...
    'Every output variable must be an N-by-1 column on tQuery.');
assert(all(cellfun(@(v) all(isfinite(v(:))), outputColumns)), ...
    'An output column contains NaN or Inf.');
assert(isscalar(maxAlphaErrorDeg) && isscalar(maxGammaErrorDeg), ...
    'Fitting-error calculations must return scalars.');
fprintf('Validated %d output columns, each with %d rows.\n', ...
    numel(outputColumns), N);

%% Save results
pitchHistory = table( ...
    tQuery, ...
    mach, ...
    altitude, ...
    velocity, ...
    dynamicPressure, ...
    rad2deg(alpha), ...
    rad2deg(gamma), ...
    rad2deg(theta), ...
    pitchRate, ...
    pitchAccel, ...
    'VariableNames', { ...
    'Time_s', ...
    'Mach', ...
    'Altitude_m', ...
    'Speed_mps', ...
    'DynamicPressure_Pa', ...
    'Alpha_deg', ...
    'FlightPath_deg', ...
    'PitchAttitude_deg', ...
    'PitchRate_rad_s', ...
    'PitchAccel_rad_s2'});

% Append aerodynamic coefficients and one clean-Cm column per candidate CG.
pitchHistory = [pitchHistory, aeroTrajectory];
for j = 1:numel(xcgCandidates)
    columnName = sprintf('CmClean_h%03d', round(100 * xcgCandidates(j)));
    pitchHistory.(columnName) = CmCleanTrajectory(:, j);
end
assert(height(pitchHistory) == N && all(isfinite(pitchHistory{:,:}), 'all'), ...
    'All final output columns must be finite and share the trajectory timeline.');

outputFile = fullfile(fileparts(mfilename('fullpath')), ...
    'pitch_maneuver_history_with_aero.csv');
writetable(pitchHistory, outputFile);
fprintf('Saved output to: %s\n', outputFile);

%% Diagnostic plots
figure;

subplot(3,1,1);
plot(tQuery, rad2deg(theta), "LineWidth", 1.5);
grid on;
ylabel("\theta (deg)");
title("Interpolated pitch motion");

subplot(3,1,2);
plot(tQuery, pitchRate, "LineWidth", 1.5);
grid on;
ylabel("q_b (rad/s)");

subplot(3,1,3);
plot(tQuery, pitchAccel, "LineWidth", 1.5);
grid on;
ylabel("dq_b/dt (rad/s^2)");
xlabel("Time (s)");


