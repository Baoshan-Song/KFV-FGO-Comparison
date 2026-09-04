%% Setup
% Reproduce the IEKF <-> Re-FGO relationship and then change one graph
% design choice at a time. Reported time is one illustrative run, not a
% benchmark. Run this file section-by-section during the tutorial.

root = fileparts(mfilename('fullpath'));
addpath(root);
environment = prepare_tutorial_environment();
root = environment.ProjectRoot;

run(fullfile(root, 'config', 'init_settings_kfv_fgo_comparison.m'));
dataPath = fullfile(root, 'data', config.data.path);
data = validate_comparison_data(load(dataPath));
outputDir = environment.OutputDirectory;

fprintf('\nKFV-FGO hands-on tutorial\n');
fprintf('Environment: MATLAB %s; local MATLAB and MATLAB Online supported.\n', ...
    environment.Release);
fprintf('Dataset: %s\n', dataPath);
fprintf('Timing note: every time value is one illustrative run, not a benchmark.\n');

%% IEKF vs Re-FGO equivalence
config.KFV.mode = 'iEKF';
config.KFV.max_iteration = 5;
config.KFV.robust_kernel = 'none';
config.KFV.window_size = 1;

fprintf('\n[1] IEKF vs Re-FGO equivalence\n');
fprintf('KFV=iEKF, max_iteration=5, kernel=none, window_size=1, imitate_KFV=1\n');

[iekfResult, iekfMetrics] = runKfvCase(config, data);
kfvEstimator = KfvEstimator(config, data);
equivalenceConfig = kfvEstimator.convert_KFV_config_to_FGO();
equivalenceConfig.FGO.autoDiff = 0;
equivalenceConfig.FGO.imitate_KFV = 1;
[refgoResult, refgoMetrics] = runFgoCase(equivalenceConfig, data);

equivalenceTable = metricsTable({'IEKF', 'Re-FGO'}, ...
    [iekfMetrics, refgoMetrics]);
disp(equivalenceTable);
maxTrajectoryDelta = max(vecnorm( ...
    iekfResult.X(1:2, :) - refgoResult.X(1:2, :), 2, 1));
fprintf('Maximum IEKF/Re-FGO horizontal trajectory difference: %.12g m\n', ...
    maxTrajectoryDelta);

equivalenceFigure = comparisonFigure(data, ...
    {iekfResult, refgoResult}, {'IEKF', 'Re-FGO'}, ...
    'IEKF and Re-FGO equivalence');
equivalenceFigureFile = fullfile(outputDir, '01_iekf_refgo_equivalence.fig');
savefig(equivalenceFigure, equivalenceFigureFile);

tutorial_results.equivalence = struct( ...
    'KFVConfiguration', config.KFV, ...
    'Configuration', equivalenceConfig, ...
    'Table', equivalenceTable, ...
    'KFVResult', iekfResult, ...
    'FGOResult', refgoResult, ...
    'MaxTrajectoryDelta', maxTrajectoryDelta);

%% Break graph memory
% Change the graph from a one-state recursive graph to a five-state
% sliding window. imitate_KFV must be disabled for a genuine SW-FGO run.
window5Iter1Config = equivalenceConfig;
window5Iter1Config.FGO.window_size = 5;
window5Iter1Config.FGO.imitate_KFV = 0;
window5Iter1Config.FGO.max_iteration = 1;

fprintf('\n[2] Graph memory: five-state SW-FGO\n');
fprintf('window_size=5, imitate_KFV=0, max_iteration=1, kernel=none\n');
[window5Iter1Result, window5Iter1Metrics] = ...
    runFgoCase(window5Iter1Config, data);
window5Iter1Table = metricsTable({'SW-FGO (w=5, iter=1)'}, ...
    window5Iter1Metrics);
disp(window5Iter1Table);

tutorial_results.window5_iter1 = struct( ...
    'Configuration', window5Iter1Config, ...
    'Table', window5Iter1Table, ...
    'Result', window5Iter1Result);

%% Change nonlinear relinearization
% Keep the five-state graph fixed and change only max_iteration.
window5Iter5Config = window5Iter1Config;
window5Iter5Config.FGO.max_iteration = 5;

fprintf('\n[3] Nonlinear relinearization\n');
fprintf('window_size=5, imitate_KFV=0, max_iteration: 1 -> 5, kernel=none\n');
[window5Iter5Result, window5Iter5Metrics] = ...
    runFgoCase(window5Iter5Config, data);
window5Iter5Table = metricsTable({'SW-FGO (w=5, iter=5)'}, ...
    window5Iter5Metrics);
disp(window5Iter5Table);

iterationFigure = comparisonFigure(data, ...
    {window5Iter1Result, window5Iter5Result}, ...
    {'SW-FGO iter=1', 'SW-FGO iter=5'}, ...
    'Five-state SW-FGO: nonlinear relinearization');
iterationFigureFile = fullfile(outputDir, '02_window5_iterations.fig');
savefig(iterationFigure, iterationFigureFile);

tutorial_results.window5_iter5 = struct( ...
    'Configuration', window5Iter5Config, ...
    'Table', window5Iter5Table, ...
    'Result', window5Iter5Result);

%% Optional robust loss
% Keep graph memory and relinearization fixed; change only the loss.
window5HuberConfig = window5Iter5Config;
window5HuberConfig.FGO.robust_kernel = 'huber';

fprintf('\n[4] Optional robust-loss challenge\n');
fprintf('window_size=5, imitate_KFV=0, max_iteration=5, kernel: none -> huber\n');
[window5HuberResult, window5HuberMetrics] = ...
    runFgoCase(window5HuberConfig, data);
window5HuberTable = metricsTable({'SW-FGO (w=5, iter=5, Huber)'}, ...
    window5HuberMetrics);
disp(window5HuberTable);

robustFigure = comparisonFigure(data, ...
    {window5Iter5Result, window5HuberResult}, ...
    {'No robust loss', 'Huber'}, ...
    'Five-state SW-FGO: robust loss');
robustFigureFile = fullfile(outputDir, '03_window5_robust_loss.fig');
savefig(robustFigure, robustFigureFile);

tutorial_results.window5_huber = struct( ...
    'Configuration', window5HuberConfig, ...
    'Table', window5HuberTable, ...
    'Result', window5HuberResult);

%% Summary
summaryTable = metricsTable( ...
    {'IEKF'; 'Re-FGO'; 'SW-FGO w=5 iter=1'; ...
    'SW-FGO w=5 iter=5'; 'SW-FGO w=5 iter=5 Huber'}, ...
    [iekfMetrics, refgoMetrics, window5Iter1Metrics, ...
    window5Iter5Metrics, window5HuberMetrics]);
disp(summaryTable);
fprintf(['Takeaway: graph memory, nonlinear relinearization, and robust ' ...
    'loss are separate design choices.\n']);
fprintf('Timing values above are illustrative single runs, not benchmarks.\n');

summaryFigure = comparisonFigure(data, ...
    {iekfResult, refgoResult, window5Iter1Result, ...
    window5Iter5Result, window5HuberResult}, ...
    {'IEKF', 'Re-FGO', 'SW-FGO iter=1', 'SW-FGO iter=5', ...
    'SW-FGO Huber'}, 'Hands-on summary');
summaryFigureFile = fullfile(outputDir, '04_hands_on_summary.fig');
savefig(summaryFigure, summaryFigureFile);

tutorial_results.summaryTable = summaryTable;
tutorial_results.figureFiles = {equivalenceFigureFile; iterationFigureFile; ...
    robustFigureFile; summaryFigureFile};
tutorial_results.environment = environment;
save(fullfile(outputDir, 'tutorial_results.mat'), 'tutorial_results');

function [result, metrics] = runKfvCase(config, data)
startTime = tic;
estimator = KfvEstimator(config, data);
result = estimator.run();
elapsed = toc(startTime);
metrics = compute_comparison_metrics(result.X, data.true_positions, elapsed);
end

function [result, metrics] = runFgoCase(config, data)
startTime = tic;
estimator = FgoEstimator(config, data);
result = estimator.run();
elapsed = toc(startTime);
metrics = compute_comparison_metrics(result.X, data.true_positions, elapsed);
end

function output = metricsTable(names, metrics)
names = string(names(:));
output = table(names, [metrics.RMSE].', [metrics.P95].', ...
    [metrics.TimeSeconds].', 'VariableNames', ...
    {'Algorithm', 'RMSE_m', 'P95_m', 'Time_s'});
end

function fig = comparisonFigure(data, results, labels, titleText)
fig = figure('Name', titleText, 'Color', 'w');
layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

trajectoryAxes = nexttile(layout, 1);
hold(trajectoryAxes, 'on');
scatter(trajectoryAxes, data.emitter_positions(1, :), ...
    data.emitter_positions(2, :), 45, 'k', 'filled', ...
    'DisplayName', 'Anchors');
plot(trajectoryAxes, data.true_positions(1, :), ...
    data.true_positions(2, :), 'k-', 'LineWidth', 1.5, ...
    'DisplayName', 'Ground truth');
colors = lines(numel(results));
for index = 1:numel(results)
    plot(trajectoryAxes, results{index}.X(1, :), results{index}.X(2, :), ...
        'LineWidth', 1.25, 'Color', colors(index, :), ...
        'DisplayName', labels{index});
end
axis(trajectoryAxes, 'equal');
grid(trajectoryAxes, 'on');
xlabel(trajectoryAxes, 'X position (m)');
ylabel(trajectoryAxes, 'Y position (m)');
title(trajectoryAxes, 'Horizontal trajectory');
legend(trajectoryAxes, 'Location', 'best');

cdfAxes = nexttile(layout, 2);
hold(cdfAxes, 'on');
for index = 1:numel(results)
    metrics = compute_comparison_metrics(results{index}.X, ...
        data.true_positions, 0);
    [probability, errorValue] = ecdf(metrics.Errors);
    plot(cdfAxes, errorValue, probability, 'LineWidth', 1.5, ...
        'Color', colors(index, :), 'DisplayName', labels{index});
end
grid(cdfAxes, 'on');
xlabel(cdfAxes, 'Horizontal error (m)');
ylabel(cdfAxes, 'Cumulative probability');
title(cdfAxes, 'Horizontal error CDF');
legend(cdfAxes, 'Location', 'best');
title(layout, titleText);
end
