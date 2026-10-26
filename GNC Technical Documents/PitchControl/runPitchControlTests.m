function results = runPitchControlTests()
%RUNPITCHCONTROLTESTS Run the numerical and export tests plus Code Analyzer.
folder = fileparts(mfilename('fullpath'));
results = runtests(fullfile(folder,'tests'));
assertSuccess(results);
files = dir(fullfile(folder,'**','*.m'));
messages = 0;
for k = 1:numel(files)
    file = fullfile(files(k).folder,files(k).name);
    issues = checkcode(file,'-id');
    for j = 1:numel(issues)
        fprintf('%s:%d [%s] %s\n',files(k).name,issues(j).line,issues(j).id,issues(j).message);
    end
    messages = messages+numel(issues);
end
fprintf('Code Analyzer: %d message(s) across %d new MATLAB files.\n',messages,numel(files));
end
