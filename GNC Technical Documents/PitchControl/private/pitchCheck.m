function row = pitchCheck(name, status, detail)
%PITCHCHECK Consistent, human-readable numerical/engineering check record.
row = table(string(name), string(status), string(detail), ...
    'VariableNames', {'Check','Status','Details'});
end
