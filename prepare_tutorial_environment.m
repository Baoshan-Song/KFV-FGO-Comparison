function environment = prepare_tutorial_environment()
%PREPARE_TUTORIAL_ENVIRONMENT Validate and initialize the hands-on tutorial.
% The setup uses paths relative to this repository so the same workflow can
% run from a local clone or from a repository opened in MATLAB Online.

projectRoot = fileparts(mfilename('fullpath'));
if verLessThan('matlab', '9.15')
    error('KFVFGO:UnsupportedMATLABRelease', ...
        'This tutorial requires MATLAB R2023b or later.');
end

requiredFiles = {
    fullfile(projectRoot, 'config', 'init_settings_kfv_fgo_comparison.m')
    fullfile(projectRoot, 'core', 'estimator', 'KfvEstimator.m')
    fullfile(projectRoot, 'core', 'estimator', 'FgoEstimator.m')
    fullfile(projectRoot, 'data', 'circle_cv_gmm_L4.mat')};
missing = requiredFiles(~cellfun(@isfile, requiredFiles));
if ~isempty(missing)
    error('KFVFGO:MissingTutorialFile', ...
        'A required tutorial file is missing: %s', missing{1});
end

addpath(genpath(projectRoot));
outputDirectory = fullfile(projectRoot, 'tutorial_output');
if ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end

probeFile = tempname(outputDirectory);
[fileId, message] = fopen(probeFile, 'w');
if fileId < 0
    error('KFVFGO:OutputNotWritable', ...
        'The tutorial output folder is not writable: %s', message);
end
fclose(fileId);
delete(probeFile);

environment = struct( ...
    'ProjectRoot', projectRoot, ...
    'OutputDirectory', outputDirectory, ...
    'Release', version('-release'), ...
    'ExecutionTarget', ...
        'MATLAB R2023b or later, including the current MATLAB Online release', ...
    'IsWritable', true);
end
