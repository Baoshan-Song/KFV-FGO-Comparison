function [X_est, P_est, X_pred, P_pred, debug_info] = miekf(X_prev, P_prev, dt, omega, f, F, Q,  toa_measurements, emitter_positions, h, H, R_base, max_iter, thres, observation_function)
if nargin < 15
    observation_function = [];
end
% Initialize debug information structure
debug_info = struct();

% Propagation
tic;  % Start timer for the prediction step
X_pred = f(X_prev, dt,omega);
P_pred = F(X_prev, dt,omega) * P_prev * F(X_prev, dt,omega)' + Q;
debug_info.prediction_time = toc;  % Record prediction time


% Measurement Update


state_size = length(X_pred);

X_est = X_pred;
P_est = P_pred;
I = eye(state_size);
X_last_est = X_est;

jacobian_all = [];
residual_norm_all = [];
tic;  % Start timer for the update step

for iter = 1:max_iter

    [y, H_all, R, jacobian_all, residual_norm_all] = evaluate_measurement_model( ...
        X_last_est, toa_measurements, emitter_positions, h, H, R_base, observation_function);
    if isempty(observation_function)
        % Preserve the legacy debug-output convention for simulations.
        jacobian_all = zeros(size(H_all));
        residual_norm_all = [];
    end
    S = H_all * P_pred * H_all' + R;
    K = P_pred * H_all' / S;                 % or: K = P_pred * H' * (S \ eye(size(S)))
    X_est = X_last_est + K * y;
    P_est = (I - K * H_all) * P_pred;    % or P_est = (I - K * H) * P_pred * (I - K * H)' + K * R * K';  % Joseph form

    % disp(norm(X_est-X_last_est)/length(X_est));
    if norm(X_est-X_last_est)/length(X_est)<thres
        break;
    end
    X_last_est = X_est;

end

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
mem_info = memory;
debug_info.memory_usage = mem_info.MemUsedMATLAB;  % Current memory usage in bytes
% Optional: Add additional metrics like CPU usage, if desired (requires external tools)
debug_info.store_time = toc;  % Record update time

end

