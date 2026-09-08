function tests = test_state_dimension_and_models
tests = functiontests(localfunctions);
end

function setupOnce(test_case)
repository_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repository_root));
test_case.TestData.repository_root = repository_root;
end

function testLegacyFourStateKfvAndFgoRemainEquivalent(test_case)
current_directory = pwd;
cleanup = onCleanup(@() cd(current_directory));
cd(test_case.TestData.repository_root);

config = struct();
run(fullfile('config', 'init_settings_kfv_fgo_comparison.m'));
data = load(fullfile('data', config.data.path));

kfv_estimator = KfvEstimator(config, data);
kfv_result = kfv_estimator.run();
fgo_config = kfv_estimator.convert_KFV_config_to_FGO();
fgo_config.FGO.autoDiff = config.FGO.autoDiff;
fgo_estimator = FgoEstimator(fgo_config, data);
fgo_result = fgo_estimator.run();

verifyEqual(test_case, config.state_dim, 4);
verifySize(test_case, kfv_result.X, [4, data.num_steps]);
verifyLessThan(test_case, max(abs(kfv_result.X - fgo_result.X), [], 'all'), 1e-9);
clear cleanup;
end

function testRealConfigurationIsAScript(test_case)
current_directory = pwd;
cleanup = onCleanup(@() cd(current_directory));
cd(test_case.TestData.repository_root);

config = struct();
run(fullfile('config', 'init_settings_gnss_ins.m'));

verifyEqual(test_case, config.data.mode, 'real');
verifyEqual(test_case, config.data.type, 'gnss_imu');
verifyEqual(test_case, config.state_dim, 10);
verifySize(test_case, config.KFV.P0, [10, 10]);
verifySize(test_case, config.KFV.Q, [10, 10]);
verifyClass(test_case, config.KFV.propagation_model, 'ImuPropagationModel');
verifyClass(test_case, config.KFV.observation_model, 'GnssObservationModel');
clear cleanup;
end

function testTenStateDatasetUsesDynamicObservation(test_case)
config.state_dim = 10;
config.KFV.mode = 'EKF';
config.KFV.dt = 1;
config.KFV.omega = [];
config.KFV.errX0 = zeros(10, 1);
config.KFV.P0 = eye(10);
config.KFV.Q = eye(10) * 1e-3;
config.KFV.f = @(x, ~, input) x + [0.01 * input; zeros(9, 1)];
config.KFV.F = @(x, dt, input) eye(length(x));
config.KFV.h = [];
config.KFV.H = [];
config.KFV.R = [];
config.KFV.observation_function = @linearObservation;
config.KFV.max_iteration = 3;
config.KFV.thres_iteration = 1e-9;
config.KFV.robust_kernel = 'none';
config.KFV.robust_delta = 1;
config.KFV.window_size = 1;
config.FGO.imitate_KFV = 1;
config.FGO.autoDiff = 0;

data = MockModelDataset();
kfv_estimator = KfvEstimator(config, data);
kfv_result = kfv_estimator.run();
fgo_config = kfv_estimator.convert_KFV_config_to_FGO();
fgo_config.FGO.autoDiff = config.FGO.autoDiff;
fgo_estimator = FgoEstimator(fgo_config, data);
fgo_result = fgo_estimator.run();

verifySize(test_case, kfv_result.X, [10, data.num_steps]);
verifyTrue(test_case, all(isfinite(kfv_result.X), 'all'));
verifyLessThan(test_case, max(abs(kfv_result.X - fgo_result.X), [], 'all'), 1e-10);
end

function testStrictImuPropagationConvention(test_case)
model = ImuPropagationModel(eye(10), 9.81);
state = [6378137; 0; 0; 1; 2; 3; zeros(3, 1); 7];
measurement.time = [100; 102];
measurement.acc = [0, 0, 9.81; 5, 5, 5];
measurement.quat = [0, 0, 0, 1; 0, 0, 0, 1];

[predicted, F] = model.computePropagation(state, measurement);

verifyEqual(test_case, predicted(1:3), state(1:3) + 2 * state(4:6), ...
    'AbsTol', 1e-10);
verifyEqual(test_case, predicted(4:6), state(4:6), 'AbsTol', 1e-10);
verifyEqual(test_case, predicted(7:10), state(7:10), 'AbsTol', 1e-12);
verifyEqual(test_case, F(1:3, 4:6), 2 * eye(3), 'AbsTol', 1e-12);
end

function [residual, H, R] = linearObservation(state, measurement)
H = measurement.H;
R = measurement.R;
residual = measurement.z - H * state;
end
