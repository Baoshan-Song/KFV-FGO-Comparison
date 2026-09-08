function [residual, H_all, R, jacobian_all, residual_norm_all] = ...
    evaluate_measurement_model(x, measurement, emitter_positions, h, H, R_config, observation_function)
%EVALUATE_MEASUREMENT_MODEL Evaluate either a model callback or legacy ranges.

if nargin >= 7 && ~isempty(observation_function)
    [residual, H_all, R] = observation_function(x, measurement);
    residual = residual(:);
    measurement_count = length(residual);
    R = measurement_covariance_matrix(R, measurement_count);

    if size(H_all, 1) ~= measurement_count || size(H_all, 2) ~= length(x)
        error('KFV:InvalidObservationJacobian', ...
            'Observation Jacobian must be %d-by-%d.', measurement_count, length(x));
    end
else
    measurement_count = size(emitter_positions, 2);
    z = measurement(:);
    H_all = zeros(measurement_count, length(x));
    predicted = zeros(measurement_count, 1);

    for i = 1:measurement_count
        H_all(i, :) = H(x, emitter_positions(:, i));
        predicted(i) = h(x, emitter_positions(:, i));
    end

    residual = z - predicted;
    R = measurement_covariance_matrix(R_config, measurement_count);
end

jacobian_all = H_all;
residual_norm_all = residual.^2;
end
