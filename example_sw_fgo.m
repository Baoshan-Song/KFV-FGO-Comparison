% TOA_CV_Demo_1
% SCRIPT Tightly coupled TOA/CV demo:
%   example (constant velocity motion in a circle (low speed))
%
% Software for use with "FGO Myth Busters".
%
% Created 05/25/2025 by Baoshan Song

% Copyright 2025, IPNL
% License: BSD; see license.txt for details


clc; clear; close all; profile off;

%% Step1: Load config
init_settings_swfgo;

%% Step 2: Generate/load data
if strcmp(config.data.mode, 'sim')
    if ~strcmp(config.data.path,'')
        data = load(config.data.path);
    else
        data = generate_data(config.sim);
    end
elseif strcmp(config.data.mode, 'real')
    data = load_real_data(config.data.path);
else
    disp('Failed to generate/load data!');
end


%% Step 3: State estimation
result1 = []; result2 = [];

% KFV
% estimator = KfvEstimator(config, data);
% result1 = estimator.run();

% FGO template
% fgo_template_config = estimator.convert_KFV_config_to_FGO();
fgo_template_config = config;
fgo_template_config.FGO.autoDiff =0;
fgo_template_config.FGO.max_iteration =1;
fgo_template_config.FGO.robust_kernel = 'none';
fgo_template_config.FGO.window_size =1;
estimator_fgo = FgoEstimator(fgo_template_config, data);
result2 = estimator_fgo.run();


%% Step 4: Results evaluation

% Statistics
% Calculate position error for KFV and FGO
% position_error_kfv = sqrt((result1.X(1,:) - data.true_positions(1,:)).^2 + ...
%     (result1.X(2,:) - data.true_positions(2,:)).^2);
position_error_fgo = sqrt((result2.X(1,:) - data.true_positions(1,:)).^2 + ...
    (result2.X(2,:) - data.true_positions(2,:)).^2);

% Mean Squared Error (MSE)
% mse_kfv = mean(position_error_kfv.^2);
mse_fgo = mean(position_error_fgo.^2);

% Root Mean Squared Error (RMSE)
% rmse_kfv = sqrt(mse_kfv);
rmse_fgo = sqrt(mse_fgo);

% Mean Absolute Error (MAE)
% mae_kfv = mean(position_error_kfv);
mae_fgo = mean(position_error_fgo);

% Maximum Error
% max_error_kfv = max(position_error_kfv);
max_error_fgo = max(position_error_fgo);

% Compute 95% absolute error (i.e. 95% errors will not surpass this)
% abs_error_95_kfv = prctile(position_error_kfv, 95);
abs_error_95_fgo = prctile(position_error_fgo, 95);

% Display Results
% disp('KFV Estimator Statistics:');
% disp(['MSE: ', num2str(mse_kfv)]);
% disp(['RMSE: ', num2str(rmse_kfv)]);
% disp(['MAE: ', num2str(mae_kfv)]);
% disp(['Max Error: ', num2str(max_error_kfv)]);
% disp(['95% Absolute Error for KFV: ', num2str(abs_error_95_kfv)]);

disp('FGO Estimator Statistics:');
disp(['MSE: ', num2str(mse_fgo)]);
disp(['RMSE: ', num2str(rmse_fgo)]);
disp(['MAE: ', num2str(mae_fgo)]);
disp(['Max Error: ', num2str(max_error_fgo)]);
disp(['95% Absolute Error for FGO: ', num2str(abs_error_95_fgo)]);


% Plot
figure; 
scatter(data.emitter_positions(1, :), data.emitter_positions(2, :), 100, 'k', 'filled'); hold on;
plot(data.true_positions(1,:)', data.true_positions(2,:)', 'b.-', 'LineWidth', 1.5, 'Marker', 'o','MarkerSize', 4); hold on;
% plot(result1.X(1,:), result1.X(2,:), 'r.-', 'LineWidth', 1.5, 'MarkerSize', 10);
plot(result2.X(1,:), result2.X(2,:), 'g.--', 'LineWidth', 1.5, 'Marker', '^', 'MarkerSize', 3);
legend('Anchor Points','True Trajectory', 'SW-FGO template traj'); % 'Noisy position Measurements',
xlabel('X Position (m)'); ylabel('Y Position (m)');
title('2D PDR Trajectory and Anchors'); axis equal;
        xlim([50,150]);
        ylim([-25, 85]);
% 
% figure;
% % plot(1:data.num_steps, position_error_kfv, 'r-', 'LineWidth', 1.5);
% % hold on;
% plot(1:data.num_steps, position_error_fgo, 'g--', 'LineWidth', 1.5);
% legend( 'FGO Position Error');
% xlabel('Time Step');
% ylabel('Position Error [m]');
% title('Position Error vs Time');
% grid on;

% figure;
% plot(1:data.num_steps, position_error_kfv-position_error_fgo, 'r-', 'LineWidth', 1.5);
% legend('KFV-FGO Position Error');
% xlabel('Time Step');
% ylabel('Position Error [m]');
% title('Position Error vs Time');
% grid on;

