% circle_toa_ekf_simulation.m
% Simulates circular trajectory + TOA measurements + EKF estimation 
% with non-Gaussian noise and nonlinearity analysis

rng(42);

%% Parameter Settings
clear; clc; close all;

% Simulation parameters
num_steps = 100;         % Number of time steps
dt = 1.0;                % Sampling interval
radius = 100;            % Radius of circular trajectory [m]
omega = 2*pi/(num_steps*dt);  % Constant angular velocity [rad/s]

% Emitter parameters
num_emitters = 6;
emitter_radius = 105;    % Radius at which emitters are placed [m]

% TOA noise parameters (Gaussian Mixture Model)
gmm_weights = [1, 0];       % Mixing weights
gmm_sigmas = [0.1, 10];     % Standard deviations [m]

%% Simulate Ground Truth Circular Trajectory
angles = linspace(0, 2*pi, num_steps);
true_positions = [radius * cos(angles); radius * sin(angles)];
true_velocities = [angles+pi/2; ones(size(angles))*omega];

%% Define Emitter Locations
emitter_angles = linspace(0, 2*pi, num_emitters+1); emitter_angles(end) = [];
emitter_positions = emitter_radius * [cos(emitter_angles); sin(emitter_angles)];

%% Simulate TOA Range Measurements
toa_measurements = zeros(num_emitters, num_steps);
for t = 1:num_steps
    for i = 1:num_emitters
        d = norm(true_positions(:,t) - emitter_positions(:,i));
        % Sample from GMM noise
        k = find(rand <= cumsum(gmm_weights),1);
        noise = gmm_sigmas(k) * randn();
        toa_measurements(i,t) = d + noise;
    end
end

% Save simulation data
save('circle_cv_gmm_L2.mat');

%% EKF Initialization
Q = diag([1e-2, 1e-2, 1e-4, 1e-4]);   % Process noise covariance
R_base = 1^2;                         % Measurement noise variance [m^2]

x_est = [200; 0; 0; omega];           % Initial state estimate
P = diag([50, 50, 0.1, 0.1]);         % Initial state covariance

% Motion model (with true angular velocity)
f = @(x, dt, omega) [
    x(1) + x(3)*dt;
    x(2) + x(4)*dt;
    x(3) - omega * x(4) * dt;
    x(4) + omega * x(3) * dt
];

% Jacobian of motion model
F = @(x, dt, omega) [
    1, 0, dt, 0;
    0, 1, 0, dt;
    0, 0, 1, -omega * dt;
    0, 0, omega * dt, 1
];

% TOA measurement model and Jacobian
h = @(x, ep) sqrt((x(1)-ep(1)).^2 + (x(2)-ep(2)).^2);
H = @(x, ep) [(x(1)-ep(1))/h(x,ep), (x(2)-ep(2))/h(x,ep), 0, 0];

x_estimates = zeros(4, num_steps);
x_estimates(:,1) = x_est;

frobenius_norms = zeros(1, num_steps);  % Nonlinearity indicator

%% EKF Loop
for t = 2:num_steps
    % True velocity (used only for reference)
    v_true = [ -radius * sin(angles(t)); radius * cos(angles(t)) ] * omega;

    % Prediction step
    x_pred = f(x_est, dt, omega);
    Fk = F(x_est, dt, omega);
    P_pred = Fk * P * Fk' + Q;

    % Update step with multiple TOA measurements
    z = toa_measurements(:,t);
    H_all = zeros(num_emitters, 4);
    h_all = zeros(num_emitters, 1);
    R = R_base * eye(num_emitters);

    for i = 1:num_emitters
        H_all(i,:) = H(x_pred, emitter_positions(:,i));
        h_all(i) = h(x_pred, emitter_positions(:,i));
    end

    y = z - h_all;
    S = H_all * P_pred * H_all' + R;
    K = P_pred * H_all' / S;
    x_est = x_pred + K * y;
    P = (eye(4) - K * H_all) * P_pred;

    x_estimates(:,t) = x_est;

    % Compute Frobenius norm of approximate Hessian
    H_f = zeros(4,4);
    for i = 1:num_emitters
        Hi = H(x_est, emitter_positions(:,i));
        H_f = H_f + Hi' * Hi;
    end
    frobenius_norms(t) = norm(H_f, 'fro');
end

%% Visualization
figure;
subplot(2,2,1);
plot(true_positions(1,:), true_positions(2,:), 'g', 'LineWidth',2); hold on;
plot(x_estimates(1,:), x_estimates(2,:), 'b--', 'LineWidth',1.5);
scatter(emitter_positions(1,:), emitter_positions(2,:), 80, 'r', 'filled');
legend('Ground Truth','EKF Estimate','Emitters'); axis equal;
title('Trajectory Comparison'); grid on;

subplot(2,2,2);
plot(1:num_steps, vecnorm(x_estimates(1:2,:) - true_positions));
title('Position Error Norm'); grid on; xlabel('Time Step'); ylabel('Error [m]');

subplot(2,2,3);
plot(2:num_steps, frobenius_norms(2:end));
title('Nonlinearity (Frobenius Norm of Hessian Approximation)');
xlabel('Time Step'); grid on;
