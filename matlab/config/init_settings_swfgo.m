%% Data configuration
% Input data source
config.data.mode = 'sim';
config.data.path = 'circle_cv_gmm_L4.mat';

%% FGO configuration: automatically tranform KFV to FGO, thus with no need for FGO configuration
% % whether using FGO template to imitate KFV (Markov assumption and closed-form solution)
config.FGO.imitate_KFV = 0;
config.FGO.autoDiff = 0;

%% FGO configuration
config.method = 'FGO';
config.FGO.dt = 1;
config.FGO.errX0 = [100,-100, 0, 0]';
config.FGO.P0 =  diag([50, 50, 1, 1]);  % initial covariance
config.FGO.radius=100;

% propagation
config.FGO.omega = 2*pi/100;
config.FGO.f = @(x, dt, omega) [
    x(1) + x(3)*dt;
    x(2) + x(4)*dt;
    x(3) - omega * x(4) * dt;
    x(4) + omega * x(3) * dt
];

config.FGO.F = @(x, dt, omega) [
    1, 0, dt, 0;
    0, 1, 0, dt;
    0, 0, 1, -omega * dt;
    0, 0, omega * dt, 1
];

config.FGO.Q = diag([1e-4, 1e-4, 1e-4, 1e-4]);   % propagation noise
% measurement 
config.FGO.h =  @(x, ep) sqrt((x(1)-ep(1)).^2 + (x(2)-ep(2)).^2);
config.FGO.H = @(x, ep) [(x(1)-ep(1))/config.KFV.h(x,ep), (x(2)-ep(2))/config.KFV.h(x,ep), 0, 0];

config.FGO.R = 0.1^2;               % m^2
% variant property
config.FGO.max_iteration = 1;
config.FGO.thres_iteration = 1E-6;
config.FGO.robust_kernel = 'huber';
config.FGO.robust_delta = 2;
config.FGO.window_size =4;



