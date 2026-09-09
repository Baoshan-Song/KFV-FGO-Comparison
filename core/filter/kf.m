function [X_est, P_est, X_pred, P_pred] = kf(X_prev, P_prev, F, Q, B, U, Z, H, R)
    % Kalman Filter step for linear system
    % X_prev: previous state estimate
    % P_prev: previous covariance
    % F: state transition
    % Q: process noise covariance
    % B: control input matrix
    % U: control input (e.g. PDR step)
    % Z: measurement
    % H: observation matrix
    % R: measurement noise covariance

    % Predict
    X_pred = F * X_prev + B * U;
    P_pred = F * P_prev * F' + Q;

    % Update
    K = P_pred * H' / (H * P_pred * H' + R);
    X_est = X_pred + K * (Z - H * X_pred);
    P_est = (eye(size(X_pred, 1)) - K * H) * P_pred;
end


