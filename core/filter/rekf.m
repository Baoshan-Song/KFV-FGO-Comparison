function [X_est, P_est, X_pred, P_pred, debug_info] = rekf(X_prev, P_prev, dt,omega,  f, F, Q, toa_measurements, emitter_positions, h, H, R_base, loss_type, delta, observation_function)
% Initialize debug information structure
debug_info = struct();

% Propagation
tic;  % Start timer for the prediction step
X_pred = f(X_prev, dt,omega);
P_pred = F(X_prev, dt,omega) * P_prev * F(X_prev, dt,omega)' + Q;
debug_info.prediction_time = toc;  % Record prediction time

tic;  % Start timer for the update step

if nargin < 14 || isempty(delta)
    delta = 1.0;
end
if nargin < 13 || isempty(loss_type)
    loss_type = 'huber';
end
if nargin < 15
    observation_function = [];
end

state_size = length(X_pred);

[y, H_all, R, jacobian_all, residual_norm_all] = evaluate_measurement_model( ...
    X_pred, toa_measurements, emitter_positions, h, H, R_base, observation_function);
S = H_all * P_pred * H_all' + R;

% Normalized residual
% r_vec = abs(y) ./ sqrt(diag(S));  % Standardized residual
r_vec = abs(y) ./ sqrt(diag(R));

% Weighting
w_vec = arrayfun(@(r) compute_weight(r, loss_type, delta), r_vec);  % m x 1
W = diag(w_vec);   % m x m

% Recompute H and y
H_r = W^(1/2) * H_all;  % m x n
y_r = W^(1/2) * y;  % m x 1
R_r = R;            % Optional, reweight H and y here is equivalent to reweighting R

% Robust kernel is only applied to residual between measurement and predict states
% Kalman Gain
S_r = H_r * P_pred * H_r' + R_r;
K = P_pred * H_r' / S_r;
X_est = X_pred + K * y_r;
P_est = (eye(state_size) - K * H_r) * P_pred;

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


function w = compute_weight(r, loss_type, delta)
switch lower(loss_type)
    case 'huber'
        if abs(r) <= delta
            w = 1;
        else
            w = delta / abs(r);
        end
    case 'cauchy'
        w = 1 / (1 + (r/delta)^2);
    case 'tukey'
        if abs(r) <= delta
            w = (1 - (r/delta)^2)^2;
        else
            w = 0;
        end
    case 'none'
        w =1;
    otherwise
        error('Unknown loss type: %s', loss_type);
end
end
