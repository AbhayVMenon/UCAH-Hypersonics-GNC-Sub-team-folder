function prepared = preparePitchTrajectory(config)
%PREPAREPITCHTRAJECTORY Smooth once and evaluate paired aero states once.
% No xCG, Iyy, or Cbar enters the interpolation. The returned structure can
% be reused by calculatePitchControlMoments(config, prepared).
config = validatePitchConfig(config);
for name = {'trajectoryFile','aeroFile'}
    key = name{1};
    assert(isfield(config,key) && isscalar(string(config.(key))) && ...
        isfile(config.(key)), 'Pitch:InputFile', 'Missing/invalid %s.', key);
end
for name = {'csaps','fnval','fnder'}
    assert(exist(name{1}, 'file') ~= 0, 'Pitch:MissingSplineToolbox', ...
        '%s is unavailable. Install/license Curve Fitting Toolbox; no fallback is used.', name{1});
end
helperFolder = fileparts(fileparts(mfilename('fullpath')));
oldPath = path;
restorePath = onCleanup(@() path(oldPath));
addpath(helperFolder);
assert(exist('interpolateAeroTrajectory','file') == 2, 'Pitch:MissingAeroHelper', ...
    'Expected existing interpolateAeroTrajectory.m in %s.', helperFolder);

raw = readtable(config.trajectoryFile, 'VariableNamingRule', 'preserve');
required = {'Time_s','Mach','Altitude_m','Speed_mps', ...
    'DynamicPressure_kPa','FlightPath_deg','Alpha_deg'};
validateColumns(raw, required, 'trajectory');
assert(height(raw) >= 3, 'Pitch:TrajectorySize', 'At least three trajectory points are required.');
sorted = sortrows(raw, 'Time_s');
t = sorted.Time_s;
assert(all(diff(t) > 0), 'Pitch:TrajectoryTime', ...
    'Trajectory time must be unique and strictly increasing after sorting.');
assert(all(sorted.DynamicPressure_kPa > 0), 'Pitch:DynamicPressure', ...
    'Trajectory DynamicPressure_kPa must be positive.');
assert(all(sorted.Mach > 0 & sorted.Speed_mps > 0), 'Pitch:TrajectoryState', ...
    'Trajectory Mach and speed must be positive.');

nIntervals = ceil((t(end)-t(1))/config.dt_s);
tq = linspace(t(1), t(end), nIntervals+1).';
alphaSpline = csaps(t, deg2rad(sorted.Alpha_deg), config.smoothParam);
gammaSpline = csaps(t, deg2rad(sorted.FlightPath_deg), config.smoothParam);
alpha = columnEval(alphaSpline,tq);
gamma = columnEval(gammaSpline,tq);
pitchRate = columnEval(fnder(alphaSpline,1),tq) + columnEval(fnder(gammaSpline,1),tq);
pitchAccel = columnEval(fnder(alphaSpline,2),tq) + columnEval(fnder(gammaSpline,2),tq);
common = table(tq, interp1(t,sorted.Mach,tq,'pchip'), ...
    interp1(t,sorted.Altitude_m,tq,'pchip'), interp1(t,sorted.Speed_mps,tq,'pchip'), ...
    1000*interp1(t,sorted.DynamicPressure_kPa,tq,'pchip'), ...
    rad2deg(alpha), rad2deg(gamma), rad2deg(alpha+gamma), pitchRate, pitchAccel, ...
    'VariableNames', {'Time_s','Mach','Altitude_m','Speed_mps','DynamicPressure_Pa', ...
    'Alpha_deg','FlightPath_deg','PitchAttitude_deg','PitchRate_rad_s','PitchAccel_rad_s2'});

aero = readtable(config.aeroFile, 'VariableNamingRule', 'preserve');
validateColumns(aero, {'Sref_m2','mass_kg'}, 'aerodynamic');
sref = constantColumn(aero.Sref_m2,'Sref_m2');
mass = constantColumn(aero.mass_kg,'mass_kg');
assert(sref > 0 && mass > 0, 'Pitch:ReferenceData', 'Database Sref and mass must be positive.');
if ~isempty(config.Sref_m2)
    assert(abs(config.Sref_m2-sref) <= 1e-10*max(abs(sref),eps), ...
        'Pitch:SrefMismatch', ['Sref override %.12g differs from dataset %.12g m^2. ' ...
        'Do not rescale moments without consistently renormalizing the coefficients.'], config.Sref_m2,sref);
end

% The existing helper validates the grid and all SMOOTHED query states.
aeroHistory = interpolateAeroTrajectory(aero,common.Mach,common.Altitude_m,common.Alpha_deg);
common = [common, aeroHistory];
nodeAero = interpolateAeroTrajectory(aero,sorted.Mach,sorted.Altitude_m,sorted.Alpha_deg);
nodes = sorted(:,required);
nodes.AlphaSmoothed_deg = rad2deg(columnEval(alphaSpline,t));
nodes.FlightPathSmoothed_deg = rad2deg(columnEval(gammaSpline,t));
nodes.AlphaFitError_deg = nodes.AlphaSmoothed_deg-sorted.Alpha_deg;
nodes.FlightPathFitError_deg = nodes.FlightPathSmoothed_deg-sorted.FlightPath_deg;
nodes.DatabaseCL = nodeAero.CL;
nodes.DatabaseCD = nodeAero.CD;
for name = {'CL','CD'}
    key = name{1};
    if ismember(key, sorted.Properties.VariableNames)
        validateColumns(sorted, {key}, 'trajectory');
        nodes.(['Trajectory' key]) = sorted.(key);
        nodes.([key 'Difference']) = nodeAero.(key)-sorted.(key);
        common.(['Trajectory' key '_pchip']) = interp1(t,sorted.(key),tq,'pchip');
    end
end
assert(all(isfinite(common{:,:}),'all'), 'Pitch:NonfiniteHistory', 'Common history contains NaN/Inf.');
assert(all(common.DynamicPressure_Pa > 0), 'Pitch:DynamicPressure', 'Interpolated pressure must be positive.');
spacingTolerance = 64*eps(max(1,max(abs(t))));
assert(max(diff(tq)) <= config.dt_s+spacingTolerance && ...
    isequal(tq([1 end]), t([1 end])), 'Pitch:TimeGrid', 'Time grid spacing/endpoints are invalid.');

axes = {'Mach','alt_m','alpha_deg'};
lo = zeros(3,1); hi = lo; counts = lo;
for k = 1:3
    values = unique(aero.(axes{k}));
    lo(k) = values(1); hi(k) = values(end); counts(k) = numel(values);
end
prepared.schemaVersion = 1;
prepared.common = common;
prepared.rawTrajectory = raw; % Original row order and all original columns.
prepared.sortedTrajectory = sorted;
prepared.nodeDiagnostics = nodes;
prepared.aeroGrid = table(string(axes(:)),counts,lo,hi, ...
    'VariableNames',{'Axis','Count','Minimum','Maximum'});
prepared.aeroSource = aero; % Enables sensitivity studies without guessing inputs.
prepared.preparationConfig = config;
prepared.sourceInfo.trajectory = dir(config.trajectoryFile);
prepared.sourceInfo.aero = dir(config.aeroFile);
prepared.metadata = struct('Sref_m2',sref,'mass_kg',mass, ...
    'maxAlphaFitError_deg',max(abs(nodes.AlphaFitError_deg)), ...
    'maxGammaFitError_deg',max(abs(nodes.FlightPathFitError_deg)), ...
    'maxTimeStep_s',max(diff(tq)), 'aeroRows',height(aero), ...
    'trajectoryRows',height(raw));
prepared.qualityChecks = [ ...
    pitchCheck('Trajectory time and required columns','PASS',sprintf('%d raw points; unique increasing sorted time.',height(raw))); ...
    pitchCheck('Common grid and endpoints','PASS',sprintf('%d points; max dt %.12g s; both endpoints retained.',height(common),max(diff(tq)))); ...
    pitchCheck('Aero grid and domain','PASS',sprintf('%d unique complete grid rows; paired lookup; no extrapolation.',height(aero))); ...
    pitchCheck('Pressure and angle units','PASS','Trajectory kPa converted once to Pa; spline derivatives calculated in radians.'); ...
    pitchCheck('Spline fit error','REVIEW',sprintf('Max alpha error %.6g deg; gamma error %.6g deg. Inspect fit and acceleration plots.',prepared.metadata.maxAlphaFitError_deg,prepared.metadata.maxGammaFitError_deg)); ...
    pitchCheck('Database references','PASS',sprintf('Sref %.12g m^2; mass %.12g kg (metadata only).',sref,mass))];
for key = ["CL","CD"]
    diffName = char(key+"Difference");
    if ismember(diffName,nodes.Properties.VariableNames)
        maxError = max(abs(nodes.(diffName)));
        status = "PASS";
        if maxError > config.coefficientComparisonTolerance, status = "REVIEW"; end
        prepared.qualityChecks = [prepared.qualityChecks; pitchCheck( ...
            key+" at original states",status,sprintf('Max |database - trajectory| %.6g; diagnostic tolerance %.6g.',maxError,config.coefficientComparisonTolerance))];
    else
        prepared.qualityChecks = [prepared.qualityChecks; pitchCheck( ...
            key+" at original states",'NOT_AVAILABLE','Original trajectory does not contain this optional coefficient.')];
    end
end
end

function out = columnEval(spline, times)
out = reshape(fnval(spline,times),[],1);
end

function validateColumns(T, names, description)
assert(all(ismember(names,T.Properties.VariableNames)), 'Pitch:MissingColumns', ...
    'Missing required %s columns: %s.',description,strjoin(setdiff(names,T.Properties.VariableNames),', '));
for k = 1:numel(names)
    v = T.(names{k});
    assert(isnumeric(v) && isreal(v) && all(isfinite(v)), 'Pitch:InvalidData', ...
        '%s column %s must contain finite real numeric values.',description,names{k});
end
end

function value = constantColumn(values, name)
assert(~isempty(values), 'Pitch:ReferenceData', 'Empty %s column.',name);
value = values(1);
assert(all(abs(values-value) <= 1e-10*max(abs(value),eps)), ...
    'Pitch:NonconstantReference', 'Database %s must be constant within relative tolerance 1e-10.',name);
end
