function [X_est, P_est, X_pred, P_pred, debug_info] = iekf(X_prev, P_prev, dt, f, F, Q,  toa_measurements, emitter_positions, h, H, R_base, max_iter, thres)
% Initialize debug information structure
debug_info = struct();

% Propagation
tic;  % Start timer for the prediction step
X_pred = f(X_prev, dt) ;
P_pred = F(X_prev, dt) * P_prev * F(X_prev, dt)' + Q;
debug_info.prediction_time = toc;  % Record prediction time

% Measurement Update
tic;  % Start timer for the update step
    
num_emitters = size(emitter_positions,2);

X_est = X_pred;
P_est = P_pred;
I = eye(length(X_est));
X_last_est = X_est;

    % Compute Jacobian and residual for all emitters
    jacobian_all = zeros(num_emitters, 4);  % Jacobian matrix for all emitters
    residual_norm_all = zeros(num_emitters, 1);  % Residual norms for each emitter
    

for iter = 1:max_iter
    z = toa_measurements;
    H_all = zeros(num_emitters, 4);
    h_all = zeros(num_emitters, 1);
    R = R_base * eye(num_emitters);
    for i = 1:num_emitters
        H_all(i,:) = H(X_est, emitter_positions(:,i));
        h_all(i) = h(X_est, emitter_positions(:,i));

        % Compute residual (difference between measured and predicted values)
        residual = z(i) - h_all(i);
        
        % Compute Jacobian (already computed)
        jacobian_all(i, :) = H_all(i, :);
        
        % Compute the residual norm (L2 norm)
        residual_norm_all(i) = residual^2;
        
    end
    y = z - h_all + H_all * (X_est - X_pred);
    S = H_all * P_pred * H_all' + R;
    K = P_pred * H_all' / S;                 % or: K = P_pred * H' * (S \ eye(size(S)))
    X_est = X_pred + K * y;
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

