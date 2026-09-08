function [X_est, P_est, X_pred, P_pred, debug_info] = ekf(X_prev, P_prev, dt, omega, f, F, Q, toa_measurements, emitter_positions, h, H, R_base, observation_function)

if nargin < 13
    observation_function = [];
end

% Initialize debug information structure
debug_info = struct();

% --- Propagation Step ---
tic;  % Start timer for the prediction step
X_pred = f(X_prev, dt,omega);
P_pred = F(X_prev, dt,omega) * P_prev * F(X_prev, dt,omega)' + Q;
debug_info.prediction_time = toc;  % Record prediction time

% --- Measurement Update Step ---
tic;  % Start timer for the update step

state_size = length(X_pred);

[y, H_all, R, jacobian_all, residual_norm_all] = evaluate_measurement_model( ...
    X_pred, toa_measurements, emitter_positions, h, H, R_base, observation_function);

% Compute the Kalman Gain
S = H_all * P_pred * H_all' + R;  % Innovation covariance
K = P_pred * H_all' / S;  % Kalman gain

% Update the state and covariance
X_est = X_pred + K * y;
P_est = (eye(state_size) - K * H_all) * P_pred;

debug_info.update_time = toc;  % Record update time

tic;
% --- Collect Debug Information ---
debug_info.jacobian_all = jacobian_all;  % All Jacobians
% debug_info.residual_all = residual_all;  % All residual for emitters
debug_info.residual_norm_all = residual_norm_all;  % Residual norms for each emitter
debug_info.Kalman_gain = K;  % Kalman gain
debug_info.innovation_covariance = S;  % Innovation covariance
debug_info.residual = y;  % Residual vector

% Measure memory usage
% mem_info = memory;
% debug_info.memory_usage = mem_info.MemUsedMATLAB;  % Current memory usage in bytes

% Optional: Add additional metrics like CPU usage, if desired (requires external tools)

debug_info.store_time = toc;  % Record update time

% --- Final Output ---
end
