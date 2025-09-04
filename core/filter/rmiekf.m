function [X_est, P_est, X_pred, P_pred, debug_info] = rmiekf(X_prev, P_prev, dt, omega, f, F, Q,  toa_measurements, emitter_positions, h, H, R_base, max_iter, thres, loss_type, delta)
% Initialize debug information structure
debug_info = struct();

% Propagation
tic;  % Start timer for the prediction step
X_pred = f(X_prev, dt,omega);
P_pred = F(X_prev, dt,omega) * P_prev * F(X_prev, dt,omega)' + Q;
debug_info.prediction_time = toc;  % Record prediction time

% Measurement Update
num_emitters = size(emitter_positions,2);
X_est = X_pred;
P_est = P_pred;
I = eye(length(X_est));
X_last_est = X_est;
W = [];

% Compute Jacobian and residual for all emitters
jacobian_all = zeros(num_emitters, 4);  % Jacobian matrix for all emitters
% residual_norm_all = zeros(num_emitters, 1);  % Residual norms for each emitter
residual_norm_all = [];
R = R_base * eye(num_emitters);
tic;  % Start timer for the update step

for iter = 1:max_iter
    z = toa_measurements;
    H_all = zeros(num_emitters, 4);
    h_all = zeros(num_emitters, 1);

    % collect the residual norm --sbs
    residual_sum = 0;
    for i = 1:num_emitters
        H_all(i,:) = H(X_last_est, emitter_positions(:,i));
        h_all(i) = h(X_last_est, emitter_positions(:,i));

        % for residual norm output, which could increase computing burden
        if 0
            % Compute residual (difference between measured and predicted values)
            residual = z(i) - h_all(i);

            % Compute Jacobian (already computed)
            jacobian_all(i, :) = H_all(i, :);

            % Compute the residual norm (L2 norm)
            residual_sum = residual_sum+  residual^2/R(i,i);
        end
    end

    % for residual norm output, which could increase computing burden
    if 0
        residual_norm_all = [residual_norm_all,[X_last_est;sqrt(residual_sum)]];
    end
    y = z - h_all ;
    % S = H_all * P_pred * H_all' + R;
    % K = P_pred * H_all' / S;                 % or: K = P_pred * H' * (S \ eye(size(S)))
    % X_est = X_last_est + K * y;
    % P_est = (I - K * H_all) * P_pred;    % or P_est = (I - K * H) * P_pred * (I - K * H)' + K * R * K';  % Joseph form


    % robust kernel part
    r_vec = abs(y) ./ sqrt(diag(R));

    % Weighting
    w_vec = arrayfun(@(r) compute_weight(r, loss_type, delta), r_vec);  % m x 1
    W = diag(w_vec);   % m x m

    % Recompute H and y
    H_r = W^(1/2) * H_all;  % m x n
    y_r = W^(1/2) * y;  % m x 1
    R_r = R;            % Optional, reweight H and y here is equivalent to reweighting R

    R= diag(1 ./ w_vec) * R;  % or R= R_base ./ w_vec

    % R_base=

    % Robust kernel is only applied to residual between measurement and predict states

    % Kalman Gain
    S_r = H_r * P_pred * H_r' + R_r;
    K = P_pred * H_r' / S_r;
    X_est = X_last_est + K * y_r;
    P_est = (I - K * H_r) * P_pred;


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
debug_info.innovation_covariance = S_r;  % Innovation covariance
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
