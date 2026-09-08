%% Real-data GNSS/IMU example
% SCRIPT Tightly coupled GNSS/IMU real-data example.

clc; clear; close all; profile off;

repository_root = fileparts(mfilename('fullpath'));
addpath(genpath(repository_root));

%% Step 1: Load config
init_settings_gnss_ins;

%% Step 2: Load real data
data = GnssImuDataset(config.data);

%% Step 3: State estimation
kfv_estimator = KfvEstimator(config, data);
kfv_result = kfv_estimator.run();

fgo_config = kfv_estimator.convert_KFV_config_to_FGO();
fgo_config.FGO.autoDiff = 0;
fgo_estimator = FgoEstimator(fgo_config, data);
fgo_result = fgo_estimator.run();

%% Step 4: Results summary
fprintf('\nProcessed %d GNSS epochs with %s.\n', data.num_steps, config.KFV.mode);
fprintf('KFV final ECEF position: %.3f %.3f %.3f m\n', ...
    kfv_result.X(1, end), kfv_result.X(2, end), kfv_result.X(3, end));
fprintf('FGO final ECEF position: %.3f %.3f %.3f m\n', ...
    fgo_result.X(1, end), fgo_result.X(2, end), fgo_result.X(3, end));
