function R = measurement_covariance_matrix(R_config, measurement_count)
%MEASUREMENT_COVARIANCE_MATRIX Expand and validate measurement covariance.
% A scalar preserves the original isotropic-noise behavior. A vector is
% interpreted as diagonal variances, and a matrix is used directly.

if isscalar(R_config)
    R = R_config * eye(measurement_count);
elseif isvector(R_config) && numel(R_config) == measurement_count
    R = diag(R_config(:));
elseif isequal(size(R_config), [measurement_count, measurement_count])
    R = R_config;
else
    error('KFV:InvalidMeasurementCovariance', ...
        'R must be scalar, length %d, or %d-by-%d.', ...
        measurement_count, measurement_count, measurement_count);
end
end
