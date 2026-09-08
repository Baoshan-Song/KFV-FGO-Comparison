%% Data configuration
% Input data directory. Change this path when the real data is stored
% outside the repository.
config_directory = fileparts(mfilename('fullpath'));
config.data.path = fullfile(fileparts(config_directory), 'data', 'urban_nav_deep');
config.data.mode = 'real';
config.data.type = 'gnss_imu';
config.data.imu_file = fullfile(config.data.path, 'xsens_imu.csv');
config.data.obs_file = fullfile(config.data.path, 'f9p_navi.obs');
config.data.nav_file = fullfile(config.data.path, 'brdm.rnx');
config.data.leap_seconds = 18;

config.method = 'KFV';
config.state_dim = 10;

%% KFV configuration
config.KFV.mode = 'EKF';
config.KFV.dt = 1;
config.KFV.omega = [];
config.KFV.errX0 = zeros(config.state_dim, 1);
config.KFV.P0 = diag([5, 5, 5, 5, 5, 5, 500, 500, 500, 100].^2);
config.KFV.Q = diag([0.3; 0.3; 0.3; 0.15; 0.15; 0.15; 0.01; 0.01; 0.01; 1e2].^2);
config.KFV.h = [];
config.KFV.H = [];
config.KFV.R = [];

config.KFV.max_iteration = 5;
config.KFV.thres_iteration = 1e-6;
config.KFV.robust_kernel = 'huber';
config.KFV.robust_delta = 1;
config.KFV.window_size = 1;

%% GNSS observation configuration
config.GNSS.minimum_snr = [];
config.GNSS.minimum_elevation_deg = [];
config.GNSS.Tref = 45;
config.GNSS.a = 30;
config.GNSS.A = 30;
config.GNSS.Fref = 10;

%% IMU propagation configuration
config.IMU.gravity_norm = 9.81;

propagation_model = ImuPropagationModel(config.KFV.Q, config.IMU.gravity_norm);
observation_model = GnssObservationModel(config.GNSS);
config.KFV.propagation_model = propagation_model;
config.KFV.observation_model = observation_model;
config.KFV.f = @(x, dt, input) propagation_model.stateTransition(x, dt, input);
config.KFV.F = @(x, dt, input) propagation_model.transitionJacobian(x, dt, input);
config.KFV.observation_function = ...
    @(x, measurement) observation_model.computeObservation(x, measurement);

config.FGO.imitate_KFV = 1;
config.FGO.autoDiff = 0;
