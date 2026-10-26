function results = runPitchControlExample()
%RUNPITCHCONTROLEXAMPLE Edit pitchControlConfig.m, then run this function.
config = pitchControlConfig();
results = calculatePitchControlMoments(config);
summary = results.summary;
disp(summary(summary.Extremum == "MaximumAbsolute" & ...
    summary.Quantity == "M_maneuver_required_Nm", ...
    {'CaseID','xCG_m','Iyy_kg_m2','SignedValue','Time_s','Mach'}));
end
