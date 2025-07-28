%% Data configuration
% Input data source
config.data.mode = 'sim';
config.data.path = 'circle_cv_gmm_L4.mat';

% Output motion, error, intermediate result files 
config.output_profile_name = 'results/TOA_CV_Demo_1_Motion.csv';
config.output_errors_name = 'results/TOA_CV_Demo_1_Errors.csv';
config.output_inter_name = 'results/TOA_CV_Demo_1_Inter.csv';

% Simulation settings
config.sim.motion = [];
config.sim.save = 0;
config.sim.save_path = 'data/TOA_CV_Demo_1_simulation.csv';

% (TODO) Real data settings

config.method = 'KFV';
%% KFV configuration
config.KFV.mode = 'riEKF';  % EKF, iEKF, rEKF, riEKF,swEKF
% predict
config.KFV.dt = 1;
config.KFV.errX0 = [100,-100, 0, 0]';
config.KFV.P0 = diag([50, 50, 1, 1]);  % initial covariance
% config.KFV.f = @(x, dt) [x(1)+x(3)*dt; x(2)+x(4)*dt; x(3); x(4)];
% config.KFV.F = @(x, dt) [1 0 dt 0; 0 1 0 dt; 0 0 1 0; 0 0 0 1];
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


% config.KFV.Q = diag([1e1, 1e1, 1e-2, 1e-2]);   % propagation noise
config.KFV.Q = diag([1e-4, 1e-4, 1e-4, 1e-4]);   % propagation noise

config.KFV.omega = 2*pi/100;
% update
config.KFV.h = @(x, ep) sqrt((x(1)-ep(1)).^2 + (x(2)-ep(2)).^2);
config.KFV.H = @(x, ep) [(x(1)-ep(1))/config.KFV.h(x,ep), (x(2)-ep(2))/config.KFV.h(x,ep), 0, 0];
config.KFV.R = 0.1^2; % m^2

% variant property
config.KFV.max_iteration =100;
config.KFV.thres_iteration = 1E-6;
config.KFV.robust_kernel = 'huber';
config.KFV.robust_delta = 2;
config.KFV.window_size = 1;

%% FGO configuration: automatically tranform KFV to FGO, thus with no need for FGO configuration
% % whether using FGO template to imitate KFV (Markov assumption and closed-form solution)
config.FGO.imitate_KFV = 1;
config.FGO.autoDiff = 0;

%% FGO configuration
% config.method = 'FGO';
config.FGO.dt = 1;
% config.FGO.X0 = [-100;   40;  90] * 1E-6;
config.FGO.errX0 = [100,-100, 0, 0]';
config.FGO.P0 =  diag([50, 50, 1, 1]);  % initial covariance
config.FGO.radius=100;
config.FGO.omega = 2*pi/100;
% config.FGO.f = @(x, dt) [x(1)+x(3)*dt; x(2)+x(4)*dt; x(3); x(4)];

% config.FGO.f = @(x, radius, dt) [
%     x(1) - radius * sin(x(3)) * x(4) * dt;
%     x(2) + radius * cos(x(3)) * x(4) * dt;
%     x(3) + x(4)*dt;
%     x(4)
% ];
config.FGO.f = @(x, dt, omega) [
    x(1) + x(3)*dt;
    x(2) + x(4)*dt;
    x(3) - omega * x(4) * dt;
    x(4) + omega * x(3) * dt
];
% config.FGO.F = @(x, dt) [1 0 dt 0; 0 1 0 dt; 0 0 1 0; 0 0 0 1];
% config.FGO.F = @(x, radius, dt) [
%     1, 0, -radius* cos(x(3)) * x(4) * dt, -radius* sin(x(3)) * dt;
%     0, 1, radius* sin(x(3)) * x(4) * dt,  radius* cos(x(3)) * dt;
%     0, 0, 1, dt;
%     0, 0, 0, 1
% ];
config.FGO.F = @(x, dt, omega) [
    1, 0, dt, 0;
    0, 1, 0, dt;
    0, 0, 1, -omega * dt;
    0, 0, omega * dt, 1
];

config.FGO.Q = diag([1e-4, 1e-4, 1e-4, 1e-4]);   % propagation noise
% update
config.FGO.h =  @(x, ep) sqrt((x(1)-ep(1)).^2 + (x(2)-ep(2)).^2);
config.FGO.H = @(x, ep) [(x(1)-ep(1))/config.KFV.h(x,ep), (x(2)-ep(2))/config.KFV.h(x,ep), 0, 0];

config.FGO.R = 0.1^2; % m^2
% variant property
config.FGO.max_iteration = 1;
config.FGO.thres_iteration = 1E-6;
config.FGO.robust_kernel = 'none';
config.FGO.robust_delta = 2;
config.FGO.window_size =20;




