function sensitivity = pitchControlSensitivity(prepared,config)
%PITCHCONTROLSENSITIVITY Fixed diagnostic variants, never an automatic retune.
% Each variant prepares ONE history for all cases, not one history per case.
config = validatePitchConfig(config);
base = evaluatePitchCases(prepared,config);
variantName = ["Baseline";"Half timestep"];
variantDt = [config.dt_s; config.dt_s/2];
variantSmooth = [config.smoothParam;config.smoothParam];
for p = unique(config.sensitivitySmoothParams(:)).'
    if p ~= config.smoothParam
        variantName(end+1,1) = "Smoothing "+string(p); %#ok<AGROW>
        variantDt(end+1,1) = config.dt_s; %#ok<AGROW>
        variantSmooth(end+1,1) = p; %#ok<AGROW>
    end
end
rows = cell(numel(variantName)*height(base.cases),1);
cursor = 0;
for v = 1:numel(variantName)
    cfg = config;
    cfg.dt_s = variantDt(v);
    cfg.smoothParam = variantSmooth(v);
    status = "OK";
    detail = "Diagnostic only; no change to baseline settings.";
    result = [];
    try
        if v == 1
            result = base;
        else
            result = evaluatePitchCases(preparePitchTrajectory(cfg),cfg);
        end
    catch problem
        status = "UNAVAILABLE";
        detail = string(problem.identifier)+": "+string(problem.message);
    end
    for k = 1:height(base.cases)
        cursor = cursor+1;
        id = base.cases.CaseID(k);
        baseline = base.caseHistory.M_maneuver_required_Nm(base.caseHistory.CaseID == id);
        basePeak = max(abs(baseline));
        row = base.cases(k,{'CaseID','xCG_m','Iyy_kg_m2'});
        row.Variant = variantName(v);
        row.dt_s = variantDt(v);
        row.SmoothParam = variantSmooth(v);
        row.Status = status;
        row.MaxAbsPitchAccel_rad_s2 = NaN;
        row.PeakAccelTime_s = NaN;
        row.PeakDynamicSigned_Nm = NaN;
        row.PeakDynamicMagnitude_Nm = NaN;
        row.PeakDynamicTime_s = NaN;
        row.PeakChange_Nm = NaN;
        row.PeakChange_percent = NaN;
        row.MaxAlphaFitError_deg = NaN;
        row.MaxGammaFitError_deg = NaN;
        if status == "OK"
            H = result.caseHistory(result.caseHistory.CaseID == id,:);
            [accelPeak,ia] = max(abs(H.PitchAccel_rad_s2));
            [peak,im] = max(abs(H.M_maneuver_required_Nm));
            row.MaxAbsPitchAccel_rad_s2 = accelPeak;
            row.PeakAccelTime_s = H.Time_s(ia);
            row.PeakDynamicSigned_Nm = H.M_maneuver_required_Nm(im);
            row.PeakDynamicMagnitude_Nm = peak;
            row.PeakDynamicTime_s = H.Time_s(im);
            row.PeakChange_Nm = peak-basePeak;
            if basePeak > 0, row.PeakChange_percent = 100*(peak-basePeak)/basePeak; end
            row.MaxAlphaFitError_deg = result.metadata.maxAlphaFitError_deg;
            row.MaxGammaFitError_deg = result.metadata.maxGammaFitError_deg;
        end
        row.Details = detail;
        rows{cursor} = row;
    end
end
sensitivity = vertcat(rows{:});
end
