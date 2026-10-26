function results = exportPitchControlResults(results)
%EXPORTPITCHCONTROLRESULTS Self-contained MATLAB exports; unique run directory.
config = results.config;
assert(isfield(config,'outputDirectory') && isscalar(string(config.outputDirectory)) && ...
    strlength(strtrim(string(config.outputDirectory))) > 0, ...
    'Pitch:OutputPath', 'Supply a nonblank outputDirectory.');
root = char(config.outputDirectory);
assert(~isfile(root),'Pitch:OutputPath','outputDirectory is an existing file: %s.',root);
if ~isfolder(root)
    [ok,message] = mkdir(root);
    assert(ok,'Pitch:OutputPath','Cannot create output directory: %s.',message);
end
stamp = char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'));
runDir = fullfile(root,['run_' stamp]);
suffix = 1;
while isfolder(runDir) || isfile(runDir)
    runDir = fullfile(root,sprintf('run_%s_%d',stamp,suffix));
    suffix = suffix+1;
end
[ok,message] = mkdir(runDir);
assert(ok,'Pitch:OutputPath','Cannot create run directory: %s.',message);
results.outputDirectory = string(runDir);
fprintf('Exporting CSVs and workbook (%d case-history rows)...\n',height(results.caseHistory));

tables = {results.common,results.caseHistory,results.summary,results.nodeDiagnostics};
names = {'common_trajectory_aero','case_history','summary','node_diagnostics'};
for k = 1:numel(tables)
    file = fullfile(runDir,[names{k} '.csv']);
    writetable(tables{k},file);
    back = readExportTable(file,tables{k});
    verifyRoundTrip(tables{k},back,names{k});
end
if ~isempty(results.sensitivity)
    writetable(results.sensitivity,fullfile(runDir,'interpolation_sensitivity.csv'));
end
results.qualityChecks = [results.qualityChecks; ...
    pitchCheck('CSV round trip','PASS','All four primary CSV table shapes, text IDs, and numeric values checked after writing.')];

workbook = fullfile(runDir,'pitch_control_results.xlsx');
inputs = inputTable(results);
[formulaExamples,formulaCells] = formulaExample(results);
sheetNames = {'Inputs_Verification','Raw_Trajectory','Common_Trajectory', ...
    'Case_History','Summary','Quality_Checks','Cases','Node_Diagnostics','Aero_Grid','Formula_Example'};
sheetTables = {inputs,results.rawTrajectory,results.common,results.caseHistory, ...
    results.summary,results.qualityChecks,results.cases,results.nodeDiagnostics,results.aeroGrid,formulaExamples};
if ~isempty(results.sensitivity)
    sheetNames{end+1} = 'Sensitivity';
    sheetTables{end+1} = results.sensitivity;
end
for k = 1:numel(sheetNames)
    assert(height(sheetTables{k}) <= 1048575 && width(sheetTables{k}) <= 16384, ...
        'Pitch:ExcelSize','Sheet %s exceeds Excel dimensions. Use fewer cases per exported study.',sheetNames{k});
    writetable(sheetTables{k},workbook,'Sheet',sheetNames{k},'UseExcel',false);
end
% Read the actual workbook independently of the CSVs before recording success.
for k = 1:numel(sheetNames)
    back = readExportTable(workbook,sheetTables{k},sheetNames{k});
    assert(isequal(size(back),size(sheetTables{k})), ...
        'Pitch:ExportShape','Workbook shape mismatch on %s.',sheetNames{k});
end
% Check numerical/time-history export values, not just workbook dimensions.
back = readExportTable(workbook,results.caseHistory,'Case_History');
verifyRoundTrip(results.caseHistory,back,'Workbook Case_History');
results.qualityChecks = [results.qualityChecks; pitchCheck('Workbook round trip','PASS', ...
    'All sheet dimensions and all case-history values checked. Formula samples contain explicit OOXML formulas and cached MATLAB values; native Excel recalculation is not tested by this runner.')];
index = strcmp(sheetNames,'Quality_Checks');
sheetTables{index} = results.qualityChecks;
writetable(results.qualityChecks,workbook,'Sheet','Quality_Checks','UseExcel',false);
formatPitchWorkbook(workbook,sheetNames,sheetTables,formulaCells);
writetable(results.qualityChecks,fullfile(runDir,'quality_checks.csv'));

results.plotFiles = strings(0,1);
if config.createPlots
    fprintf('Generating comparison and interpolation diagnostic plots...\n');
    results.plotFiles = plotPitchControlResults(results,runDir);
end
results.files = struct('workbook',string(workbook), ...
    'mat',string(fullfile(runDir,'pitch_control_results.mat')), ...
    'summary',string(fullfile(runDir,'summary.csv')), ...
    'history',string(fullfile(runDir,'case_history.csv')));
writeRunReadme(results);
save(results.files.mat,'results','-v7.3');
end

function verifyRoundTrip(original,back,description)
assert(isequal(size(original),size(back)), 'Pitch:ExportShape','Shape mismatch: %s.',description);
for j = 1:width(original)
    name = original.Properties.VariableNames{j};
    a = original.(name); b = back.(name);
    if isnumeric(a) || islogical(a)
        assert(isnumeric(b) || islogical(b), 'Pitch:ExportValues','Lost numeric type in %s.%s.',description,name);
        assert(all((isnan(a) & isnan(b)) | abs(double(a)-double(b)) <= 1e-10*max(1,abs(double(a)))), ...
            'Pitch:ExportValues','Numeric mismatch: %s.%s.',description,name);
    else
        assert(isequal(string(a),string(b)), 'Pitch:ExportValues','Text mismatch: %s.%s.',description,name);
    end
end
end

function T = readExportTable(file,prototype,sheet)
% Avoid type inference turning dimensionless unit text "1" into missing data.
if nargin < 3
    opts = detectImportOptions(file,'VariableNamingRule','preserve');
    opts.VariableNamesLine = 1;
    opts.VariableUnitsLine = 0;
    opts.VariableDescriptionsLine = 0;
    opts.DataLines = [2 Inf];
else
    opts = detectImportOptions(file,'Sheet',sheet,'VariableNamingRule','preserve');
end
for j = 1:width(prototype)
    name = prototype.Properties.VariableNames{j};
    v = prototype.(name);
    if isnumeric(v) || islogical(v)
        opts = setvartype(opts,name,'double');
    else
        opts = setvartype(opts,name,'string');
        opts = setvaropts(opts,name,'WhitespaceRule','preserve');
    end
end
if nargin < 3
    T = readtable(file,opts);
else
    T = readtable(file,opts,'UseExcel',false);
end
end

function T = inputTable(r)
c = r.config;
keys = ["Status";"Sref";"Cbar";"MomentReference";"mass_kg";"dt_s";"smoothParam"; ...
    "Iyy";"CG";"CoordinateDirection";"PitchSign";"trajectoryFile";"aeroFile";"CaseMode"; ...
    "Cbar discrepancy";"Model scope";"Common history";"M_vehicle_Nm";"M_inertial_Nm"; ...
    "M_trim_required_Nm";"M_maneuver_required_Nm";"Cm_CG";"DeltaCm_trim";"DeltaCm_maneuver"; ...
    "Workbook behavior";"Formula_Example";"Sign extrema";"Source mass";"Absent sensitivity values"];
numbers = nan(numel(keys),1);
numbers(2:7) = [r.metadata.Sref_m2;c.Cbar_m;c.xMomentReference_m;r.metadata.mass_kg;c.dt_s;c.smoothParam];
textValues = strings(numel(keys),1);
textValues(1) = r.status;
textValues(8) = "Per-case positive inertia; see Cases sheet";
textValues(9) = "Dimensional CG in metres; see Cases sheet";
textValues(10) = "Aft-positive x; longitudinal shift only; zero vertical offset";
textValues(11) = "Positive/negative pitch; nose-up/down mapping requires confirmation";
textValues(12:14) = [string(c.trajectoryFile);string(c.aeroFile);c.caseMode];
textValues(15) = "Handoff 0.18038 m vs legacy script 0.18188 m; new CSV does not establish Cbar";
textValues(16) = "Planar theta=alpha+gamma; constant Iyy per case; no damping/thrust or unprovided moments";
textValues(17) = "Trajectory pressure in Pa; smoothed alpha for paired linear aero lookup; no extrapolation";
textValues(18:24) = ["q*Sref*Cbar*Cm_CG";"Iyy*pitchAccel";"-M_vehicle_Nm"; ...
    "M_inertial_Nm-M_vehicle_Nm (already includes trim)"; ...
    "Cm_nose+CN*(xCG-xMomentReference)/Cbar";"-Cm_CG";"M_maneuver_required_Nm/(q*Sref*Cbar)"];
textValues(25) = "Main tables are frozen MATLAB numerical exports. Edit config and rerun for a new study; editing this workbook does not update those tables.";
textValues(26) = "Three editable sample calculations with real Excel formulas, cached MATLAB values, and automatic recalculation on open. Samples do not drive the summary.";
textValues(27) = "HasPositive/HasNegative indicate whether either sign occurs. Ties choose the first chronological row.";
textValues(28) = "Database metadata only; does not determine Iyy or establish the mass used by the trajectory optimizer.";
textValues(29) = "Blank sensitivity numbers mean unavailable variant or undefined percentage when baseline peak is zero; inspect Status and Details.";
units = strings(numel(keys),1);
units(2:9) = ["m^2";"m";"m";"kg";"s";"1";"kg m^2";"m"];
units(18:21) = "N m";
units(22:24) = "1";
sources = strings(numel(keys),1);
confirmation = repmat("NOT_APPLICABLE",numel(keys),1);
for k = 1:numel(keys)
    key = char(keys(k));
    if isfield(c.verification,key)
        sources(k) = c.verification.(key).source;
        confirmation(k) = "UNCONFIRMED";
        if c.verification.(key).confirmed, confirmation(k) = "USER_CONFIRMED"; end
    end
end
sources(5) = string(c.aeroFile)+": mass_kg";
sources(12) = string(c.trajectoryFile);
sources(13) = string(c.aeroFile);
T = table(keys,numbers,textValues,units,sources,confirmation, ...
    'VariableNames',{'Item','NumericValue','TextValue','Unit','Source','Confirmation'});
end

function [T,formulaCells] = formulaExample(r)
H = r.caseHistory(r.caseHistory.CaseID == r.cases.CaseID(1),:);
indices = unique([1;ceil(height(H)/2);height(H)]);
H = H(indices,:);
n = height(H);
c = r.config;
T = table(H.Time_s,H.xCG_m,H.Iyy_kg_m2,repmat(c.xMomentReference_m,n,1), ...
    repmat(r.metadata.Sref_m2,n,1),repmat(c.Cbar_m,n,1),H.DynamicPressure_Pa, ...
    H.PitchAccel_rad_s2,H.CN,H.Cm_nose,H.Cm_CG,H.M_vehicle_Nm,H.M_inertial_Nm, ...
    H.M_trim_required_Nm,H.M_maneuver_required_Nm,H.DeltaCm_trim,H.DeltaCm_maneuver, ...
    'VariableNames',{'Time_s','xCG_m','Iyy_kg_m2','xReference_m','Sref_m2','Cbar_m', ...
    'DynamicPressure_Pa','PitchAccel_rad_s2','CN','Cm_nose','Cm_CG','M_vehicle_Nm', ...
    'M_inertial_Nm','M_trim_required_Nm','M_maneuver_required_Nm','DeltaCm_trim','DeltaCm_maneuver'});
formulaCells = cell(n*7,2);
for k = 1:n
    row = k+1;
    expressions = {sprintf('J%d+I%d*(B%d-D%d)/F%d',row,row,row,row,row), ...
        sprintf('G%d*E%d*F%d*K%d',row,row,row,row), sprintf('C%d*H%d',row,row), ...
        sprintf('-L%d',row),sprintf('M%d-L%d',row,row),sprintf('-K%d',row), ...
        sprintf('O%d/(G%d*E%d*F%d)',row,row,row,row)};
    for j = 1:7
        index = (k-1)*7+j;
        formulaCells(index,:) = {sprintf('%c%d',char('K'+j-1),row),expressions{j}};
    end
end
end

function writeRunReadme(r)
file = fullfile(r.outputDirectory,'README.md');
[fid,message] = fopen(file,'wt');
assert(fid >= 0,'Pitch:OutputPath','Cannot write README: %s.',message);
closer = onCleanup(@() fclose(fid));
fprintf(fid,'# Pitch-control study\n\n%s\n\n',r.status);
fprintf(fid,'Generated by MATLAB %s. %d common trajectory rows, %d cases.\n\n',version,height(r.common),height(r.cases));
fprintf(fid,'Sref = %.12g m^2; Cbar = %.12g m; mass metadata = %.12g kg.\n\n',r.metadata.Sref_m2,r.config.Cbar_m,r.metadata.mass_kg);
fprintf(fid,'Sources: `%s` and `%s`.\n\n',r.config.trajectoryFile,r.config.aeroFile);
fprintf(fid,'Main workbook tables and CSVs are numerical MATLAB snapshots. Change `pitchControlConfig.m` and rerun `runPitchControlExample` for a new study. The Formula_Example sheet alone contains editable Excel formulas with cached values; it does not update the main histories or summary.\n\n');
fprintf(fid,'Moment signs mean positive/negative pitch, pending confirmation. Cbar remains disputed: handoff 0.18038 m versus old script 0.18188 m. The new CSV does not confirm Cbar or the moment reference.\n\n');
fprintf(fid,'Inspect `node_diagnostics.csv`, `quality_checks.csv`, and `interpolation_sensitivity.csv` when present. Blank sensitivity values indicate an unavailable variant or undefined relative change; inspect Status and Details. Denser output sampling does not recover dynamics missing from the original collocation points.\n\n');
fprintf(fid,'Vehicle moments are clean aerodynamic moments along the prescribed trajectory. Control requirements are not available authority, deflections, actuator loads, or proof of trim feasibility. See the source-folder README for APIs, dependencies, equations, tests, and examples.\n');
end
