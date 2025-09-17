%% Data configuration
% Input data source
config.data.mode = 'sim';
config.data.path = 'circle_cv_gmm_L4.mat';

config.method = 'KFV';
%% KFV configuration
config.KFV.mode = 'EKF';  % EKF, iEKF, rEKF, riEKF
% predict
config.KFV.dt = 1;
config.KFV.errX0 = [100,-100, 0, 0]';
% config.KFV.errX0 = [0,0, 0, 0]';
config.KFV.P0 = diag([50, 50, 1, 1]);  % initial covariance

config.KFV.f = @(x, dt, omega) [
    x(1) + x(3)*dt;
    x(2) + x(4)*dt;
    x(3) - omega * x(4) * dt;
    x(4) + omega * x(3) * dt
];
config.KFV.F = @(x, dt, omega) [
    1, 0, dt, 0;
    0, 1, 0, dt;
    0, 0, 1, -omega * dt;
    0, 0, omega * dt, 1
];

config.KFV.Q = diag([1e-4, 1e-4, 1e-4, 1e-4]);   % propagation noise

config.KFV.omega = 2*pi/100;
% update
config.KFV.h = @(x, ep) sqrt((x(1)-ep(1)).^2 + (x(2)-ep(2)).^2);
config.KFV.H = @(x, ep) [(x(1)-ep(1))/config.KFV.h(x,ep), (x(2)-ep(2))/config.KFV.h(x,ep), 0, 0];
config.KFV.R = 0.1^2; % m^2

% variant property
config.KFV.max_iteration =2;
config.KFV.thres_iteration = 1E-6;
config.KFV.robust_kernel = 'huber';
config.KFV.robust_delta = 2;
config.KFV.window_size = 1;

%% FGO configuration: automatically tranform KFV to FGO, thus with no need for FGO configuration
% % whether using FGO template to imitate KFV (Markov assumption and closed-form solution)
config.FGO.imitate_KFV = 1;
config.FGO.autoDiff = 0;

