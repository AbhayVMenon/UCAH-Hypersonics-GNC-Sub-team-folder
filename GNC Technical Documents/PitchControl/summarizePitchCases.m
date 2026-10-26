function summary = summarizePitchCases(history)
%SUMMARIZEPITCHCASES Independent signed and absolute extrema at actual rows.
quantities = ["M_vehicle_Nm","M_trim_required_Nm","M_maneuver_required_Nm", ...
    "Cm_CG","DeltaCm_trim","DeltaCm_maneuver"];
cases = unique(history.CaseID,'stable');
rows = cell(numel(cases)*numel(quantities)*3,1);
cursor = 0;
for k = 1:numel(cases)
    globalIndices = find(history.CaseID == cases(k));
    H = history(globalIndices,:);
    for quantity = quantities
        values = H.(quantity);
        [~,minRow] = min(values);
        [~,maxRow] = max(values);
        [~,absRow] = max(abs(values));
        selected = [minRow maxRow absRow];
        types = ["Minimum","Maximum","MaximumAbsolute"];
        unit = "1";
        if endsWith(quantity,"_Nm"), unit = "N m"; end
        for j = 1:3
            index = selected(j);
            cursor = cursor+1;
            row = H(index,{'CaseID','CaseLabel','xCG_m','Iyy_kg_m2'});
            row.Quantity = quantity;
            row.Extremum = types(j);
            row.Unit = unit;
            row.SignedValue = values(index);
            row.Magnitude = abs(values(index));
            row.TimeRowIndex = index;
            row.HistoryRowIndex = globalIndices(index);
            row.HasPositive = any(values > 0);
            row.HasNegative = any(values < 0);
            rows{cursor} = [row,H(index,{'Time_s','Mach','Altitude_m','Alpha_deg','DynamicPressure_Pa'})];
        end
    end
end
summary = vertcat(rows{:});
end
