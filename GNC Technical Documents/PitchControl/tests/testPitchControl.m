function tests = testPitchControl
%TESTPITCHCONTROL Physical identities, data contracts, interpolation, exports.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
folder = fileparts(fileparts(mfilename('fullpath')));
testCase.TestData.oldPath = path;
addpath(folder,fileparts(folder));
config = pitchControlConfig();
config.exportResults = false;
config.createPlots = false;
config.runSensitivity = false;
testCase.TestData.config = config;
testCase.TestData.prepared = preparePitchTrajectory(config);
end

function teardownOnce(testCase)
path(testCase.TestData.oldPath);
end

function testCommonTrajectoryAndSourceUnits(testCase)
p = testCase.TestData.prepared;
c = testCase.TestData.config;
t = p.sortedTrajectory.Time_s;
verifyEqual(testCase,p.common.Time_s([1 end]),t([1 end]));
verifyLessThanOrEqual(testCase,max(diff(p.common.Time_s)),c.dt_s+1e-12);
verifyEqual(testCase,height(p.common),ceil((t(end)-t(1))/c.dt_s)+1);
expected = 1000*interp1(t,p.sortedTrajectory.DynamicPressure_kPa,p.common.Time_s,'pchip');
verifyEqual(testCase,p.common.DynamicPressure_Pa,expected,'AbsTol',1e-10);
verifyEqual(testCase,p.metadata.Sref_m2,0.02555495,'AbsTol',1e-14);
verifyEqual(testCase,p.metadata.mass_kg,120);
verifyLessThan(testCase,max(abs(p.nodeDiagnostics.CLDifference)),1e-10);
verifyLessThan(testCase,max(abs(p.nodeDiagnostics.CDDifference)),1e-10);
end

function testReferenceCGRecoversSourceMoment(testCase)
c = testCase.TestData.config; p = testCase.TestData.prepared;
c.xCG_m = c.xMomentReference_m;
r = evaluatePitchCases(p,c);
verifyEqual(testCase,height(r.cases),1);
verifyEqual(testCase,r.caseHistory.Cm_CG,p.common.Cm_nose,'AbsTol',1e-12);
end

function testZeroAccelerationEqualsTrim(testCase)
c = testCase.TestData.config; p = testCase.TestData.prepared;
p.common.PitchAccel_rad_s2(:) = 0;
r = evaluatePitchCases(p,c);
verifyEqual(testCase,r.caseHistory.M_maneuver_required_Nm,r.caseHistory.M_trim_required_Nm);
verifyEqual(testCase,r.caseHistory.DeltaCm_maneuver,r.caseHistory.DeltaCm_trim,'AbsTol',1e-12);
end

function testInertiaSensitivityAndStaticIndependence(testCase)
c = testCase.TestData.config; p = testCase.TestData.prepared;
c.xCG_m = [1.12 1.12]; c.Iyy_kg_m2 = [10 25];
r = evaluatePitchCases(p,c);
n = height(p.common);
a = r.caseHistory(1:n,:); b = r.caseHistory(n+1:end,:);
verifyEqual(testCase,a.M_vehicle_Nm,b.M_vehicle_Nm);
verifyEqual(testCase,a.M_trim_required_Nm,b.M_trim_required_Nm);
verifyEqual(testCase,b.M_maneuver_required_Nm-a.M_maneuver_required_Nm, ...
    15*p.common.PitchAccel_rad_s2,'AbsTol',1e-9);
end

function testDimensionalCGShift(testCase)
c = testCase.TestData.config; p = testCase.TestData.prepared;
c.xCG_m = [1.08 1.13];
r = evaluatePitchCases(p,c); n = height(p.common);
a = r.caseHistory(1:n,:); b = r.caseHistory(n+1:end,:);
verifyEqual(testCase,b.Cm_CG-a.Cm_CG,p.common.CN*0.05/c.Cbar_m,'AbsTol',1e-12);
verifyEqual(testCase,b.M_vehicle_Nm-a.M_vehicle_Nm, ...
    p.common.DynamicPressure_Pa*p.metadata.Sref_m2.*p.common.CN*0.05,'AbsTol',1e-9);
end

function testScalarAndRepeatedInertiaMatch(testCase)
c = testCase.TestData.config; p = testCase.TestData.prepared;
a = evaluatePitchCases(p,c);
c.Iyy_kg_m2 = repmat(c.Iyy_kg_m2,size(c.xCG_m));
b = evaluatePitchCases(p,c);
verifyEqual(testCase,a.caseHistory,b.caseHistory);
end

function testCartesianCasesAndSharedHistory(testCase)
c = testCase.TestData.config; p = testCase.TestData.prepared;
c.xCG_m = [1.111111111111 1.222222222222]; c.Iyy_kg_m2 = [10 20 30];
c.caseMode = "cartesian";
r = evaluatePitchCases(p,c);
verifyEqual(testCase,height(r.cases),6);
verifyEqual(testCase,r.cases.xCG_m,repmat(c.xCG_m(:),3,1));
verifyEqual(testCase,r.cases.Iyy_kg_m2,repelem(c.Iyy_kg_m2(:),2));
verifyEqual(testCase,numel(unique(r.cases.CaseID)),6);
for k = 1:6
    H = r.caseHistory(r.caseHistory.CaseID == r.cases.CaseID(k),:);
    verifyEqual(testCase,H(:,p.common.Properties.VariableNames),p.common);
end
end

function testInvalidPairedDimensions(testCase)
c = testCase.TestData.config;
c.xCG_m = [1 2 3]; c.Iyy_kg_m2 = [10 20];
verifyError(testCase,@() evaluatePitchCases(testCase.TestData.prepared,c),'Pitch:PairedDimensions');
end

function testGridNodesAndArbitraryRowOrder(testCase)
aero = testCase.TestData.prepared.aeroSource;
shuffled = aero(end:-1:1,:);
r = interpolateAeroTrajectory(shuffled,aero.Mach,aero.alt_m,aero.alpha_deg);
names = {'CN','CA','CL','CD','CDpres','CDfric','Cm_nose'};
verifyEqual(testCase,r{:,names},aero{:,names},'AbsTol',1e-12);
end

function testPairedLinearInterpolationManufacturedField(testCase)
[m,h,a] = ndgrid([3 8],[0 36000],[-10 25]);
field = 2*m+0.001*h-0.3*a+0.00001*m.*h;
A = table(m(:),h(:),a(:),field(:),field(:),field(:),ones(numel(m),1), ...
    ones(numel(m),1),ones(numel(m),1),field(:),'VariableNames', ...
    {'Mach','alt_m','alpha_deg','CN','CA','CL','CD','CDpres','CDfric','Cm_nose'});
mq = [3.5 5 7.5]; hq = [2000 18000 30000]; aq = [-5 0 20];
R = interpolateAeroTrajectory(A,mq,hq,aq);
expected = 2*mq+0.001*hq-0.3*aq+0.00001*mq.*hq;
verifyEqual(testCase,height(R),3);
verifyEqual(testCase,R.CN,expected(:),'AbsTol',1e-12);
end

function testOutOfBoundsAndInvalidGrids(testCase)
A = testCase.TestData.prepared.aeroSource;
verifyError(testCase,@() interpolateAeroTrajectory(A,9,10000,0),'Aero:OutOfRange');
verifyError(testCase,@() interpolateAeroTrajectory(A(1:end-1,:),5,10000,0),'Aero:IncompleteGrid');
verifyError(testCase,@() interpolateAeroTrajectory([A;A(1,:)],5,10000,0),'Aero:DuplicateState');
end

function testAnalyticalMomentsAndIndependentExtrema(testCase)
[p,c] = smallHistory(testCase);
r = evaluatePitchCases(p,c);
verifyEqual(testCase,r.caseHistory.M_vehicle_Nm,[-4;4;8]);
verifyEqual(testCase,r.caseHistory.M_inertial_Nm,[0;6;-3]);
verifyEqual(testCase,r.caseHistory.M_trim_required_Nm,[4;-4;-8]);
verifyEqual(testCase,r.caseHistory.M_maneuver_required_Nm,[4;2;-11]);
S = r.summary;
vehicle = S(S.Quantity == "M_vehicle_Nm" & S.Extremum == "MaximumAbsolute",:);
coefficient = S(S.Quantity == "Cm_CG" & S.Extremum == "MaximumAbsolute",:);
verifyEqual(testCase,vehicle.TimeRowIndex,3);
verifyEqual(testCase,coefficient.TimeRowIndex,1);
verifyEqual(testCase,coefficient.SignedValue,-2);
for k = 1:height(S)
    H = r.caseHistory(S.HistoryRowIndex(k),:);
    verifyEqual(testCase,H.(S.Quantity(k)),S.SignedValue(k));
    verifyEqual(testCase,H.Time_s,S.Time_s(k));
    verifyEqual(testCase,abs(S.SignedValue(k)),S.Magnitude(k));
end
end

function testAbsentSignFlags(testCase)
[p,c] = smallHistory(testCase);
p.common.Cm_nose(:) = 1;
r = evaluatePitchCases(p,c);
S = r.summary(r.summary.Quantity == "M_trim_required_Nm",:);
verifyFalse(testCase,any(S.HasPositive));
verifyTrue(testCase,all(S.HasNegative));
end

function testReferenceMismatchFails(testCase)
c = testCase.TestData.config;
c.Sref_m2 = 2*testCase.TestData.prepared.metadata.Sref_m2;
verifyError(testCase,@() evaluatePitchCases(testCase.TestData.prepared,c),'Pitch:SrefMismatch');
end

function testChangedPreparationSettingsFail(testCase)
c = testCase.TestData.config;
c.dt_s = 0.1;
verifyError(testCase,@() calculatePitchControlMoments(c,testCase.TestData.prepared),'Pitch:StalePreparation');
end

function testDuplicateTrajectoryTimeFails(testCase)
c = testCase.TestData.config;
folder = tempname; mkdir(folder);
cleanup = onCleanup(@() rmdir(folder,'s'));
T = testCase.TestData.prepared.rawTrajectory;
T.Time_s(2) = T.Time_s(1);
c.trajectoryFile = fullfile(folder,'duplicate_times.csv');
writetable(T,c.trajectoryFile);
verifyError(testCase,@() preparePitchTrajectory(c),'Pitch:TrajectoryTime');
end

function testUnsortedInputPreservedAndTimeSorted(testCase)
c = testCase.TestData.config;
folder = tempname; mkdir(folder);
cleanup = onCleanup(@() rmdir(folder,'s'));
T = testCase.TestData.prepared.rawTrajectory(end:-1:1,:);
c.trajectoryFile = fullfile(folder,'reverse.csv');
writetable(T,c.trajectoryFile);
p = preparePitchTrajectory(c);
verifyEqual(testCase,p.rawTrajectory.Time_s,T.Time_s,'AbsTol',1e-10);
verifyEqual(testCase,p.common.Time_s,testCase.TestData.prepared.common.Time_s,'AbsTol',1e-10);
end

function testFormulaWorkbookAndUniqueExportFolders(testCase)
[p,c] = smallHistory(testCase);
folder = tempname; mkdir(folder);
cleanup = onCleanup(@() rmdir(folder,'s'));
c.outputDirectory = folder; c.createPlots = false;
r = evaluatePitchCases(p,c);
one = exportPitchControlResults(r);
two = exportPitchControlResults(r);
verifyNotEqual(testCase,one.outputDirectory,two.outputDirectory);
verifyTrue(testCase,isfile(one.files.mat));
verifyTrue(testCase,isfile(one.files.workbook));
T = readtable(one.files.workbook,'Sheet','Case_History','UseExcel',false);
verifyEqual(testCase,T.M_maneuver_required_Nm,[4;2;-11]);
extracted = fullfile(folder,'xml'); unzip(one.files.workbook,extracted);
text = fileread(fullfile(extracted,'xl','worksheets','sheet10.xml'));
verifyEqual(testCase,numel(regexp(text,'<f>','match')),21);
verifyTrue(testCase,contains(text,'J2+I2*(B2-D2)/F2'));
verifyTrue(testCase,contains(text,'ySplit="1"'));
end

function [p,c] = smallHistory(testCase)
p = testCase.TestData.prepared; c = testCase.TestData.config;
p.common = p.common(1:3,:);
p.common.Time_s = [0;1;2];
p.common.CN = ones(3,1);
p.common.Cm_nose = [-2;1;0.5];
p.common.DynamicPressure_Pa = [2;4;16];
p.common.PitchAccel_rad_s2 = [0;2;-1];
p.metadata.Sref_m2 = 1;
c.Cbar_m = 1; c.xMomentReference_m = 0;
c.xCG_m = 0; c.Iyy_kg_m2 = 3;
end
