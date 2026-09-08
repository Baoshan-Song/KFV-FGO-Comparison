classdef GnssPseudorangeFactor < factor
    methods
        function obj = GnssPseudorangeFactor(states, measurement, config)
            obj@factor(states, config);
            obj.z = measurement;
        end

        function obj = evaluate(obj)
            state_value = obj.states(1).value;
            [residual, H, R] = obj.config.FGO.observation_function(state_value, obj.z);
            residual = residual(:);

            standardized_residual = abs(residual) ./ sqrt(diag(R));
            weights = arrayfun(@(value) robustWeight(value, ...
                obj.config.FGO.robust_kernel, obj.config.FGO.robust_delta), ...
                standardized_residual);
            square_root_weight = diag(sqrt(weights));

            square_root_covariance = chol(R, 'lower');
            obj.A = square_root_weight * (square_root_covariance \ H);
            obj.b = square_root_weight * (square_root_covariance \ residual);
            obj.Omega = diag(weights) / R;
        end
    end
end

function weight = robustWeight(residual, loss_type, delta)
switch lower(loss_type)
    case 'huber'
        if abs(residual) <= delta
            weight = 1;
        else
            weight = delta / abs(residual);
        end
    case 'cauchy'
        weight = 1 / (1 + (residual / delta)^2);
    case 'tukey'
        if abs(residual) <= delta
            weight = (1 - (residual / delta)^2)^2;
        else
            weight = 0;
        end
    case 'none'
        weight = 1;
    otherwise
        error('KFV:UnknownRobustKernel', 'Unknown loss type: %s', loss_type);
end
end
