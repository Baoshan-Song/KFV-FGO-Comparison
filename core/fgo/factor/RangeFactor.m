classdef RangeFactor < factor
    properties

    end

    methods
        function obj = RangeFactor(states, z, Omega)
            obj@factor(states,z,Omega);
        end

        function obj = evaluate(obj)

            if obj.z.autoDiff == 0
                % Extract states, emitter position and range measurement
                state = obj.states(1);
                emitter = obj.z.emitter;
                range = obj.z.range;

                % Current geometric range
                geometric_range =  sqrt((state.value(1)-emitter(1)).^2 + (state.value(2)-emitter(2)).^2);

                residual = range - geometric_range;
                r_vec = abs(residual) .* sqrt(diag(obj.Omega));

                J =  [(state.value(1)-emitter(1))/geometric_range, ...
                    (state.value(2)-emitter(2))/geometric_range, 0, 0];

                % Apply robust kernel
                loss_type = obj.z.loss_type;  % e.g., 'huber', 'tukey'
                delta = obj.z.loss_delta;          % threshold parameter
                [~, sqrt_rho_prime] = compute_robust_weight(r_vec, loss_type, delta);

                % (TODO): mutiply the obj.Omega to A
                L = chol(obj.Omega, 'lower');
                obj.A = L' *sqrt_rho_prime* J;
                obj.b = L' *sqrt_rho_prime* residual;

                obj.Omega = obj.Omega*sqrt_rho_prime^2; % test --sbs

            elseif  obj.z.autoDiff == 1
                % Extract states, emitter position and range measurement
                state = obj.states(1);
                emitter = obj.z.emitter;
                range = obj.z.range;

                % Forward with autodiff
                x = dlarray(state.value);  % [x y vx vy]

                % Wrap cost function
                residual_fun = @(x) range - sqrt((x(1)-emitter(1))^2 + (x(2)-emitter(2))^2);

                % Compute residual and Jacobian via automatic differentiation
                [residual, grad] = dlfeval(@(x) valueAndGrad(residual_fun, x), x);

                % Convert to numeric
                residual = extractdata(residual);
                J = extractdata(grad);

                % test --sbs
                J = -J';

                % Robust kernel
                loss_type = obj.z.loss_type;
                delta = obj.z.loss_delta;

                % Robust residual (scalar residual, apply to 1x1 loss)
                r_vec = abs(residual) * sqrt(diag(obj.Omega));
                [~, sqrt_rho_prime] = compute_robust_weight(r_vec, loss_type, delta);

                % Final info matrix
                L = chol(obj.Omega, 'lower');

                % Assign to factor
                obj.A = L' * sqrt_rho_prime * J;        % [1x4]
                obj.b = L' * sqrt_rho_prime * residual; % [1x1]
                
                obj.Omega = obj.Omega*sqrt_rho_prime^2; % test --sbs

            else
                disp("Please define autoDiff!");
            end

        end


    end
end

function [y, grad] = valueAndGrad(fun, x)
    y = fun(x);
    grad = dlgradient(y, x);
end

function [w, sqrt_rho_prime] = compute_robust_weight(r, loss_type, delta)
    s = r^2;  % squared residual

    switch lower(loss_type)
        case 'huber'
            if abs(r) <= delta
                sqrt_rho_prime = 1;
            else
                sqrt_rho_prime = sqrt(delta / abs(r));  % sqrt of derivative of rho(s)
            end
        case 'cauchy'
            sqrt_rho_prime = sqrt(1 / (1 + (r/delta)^2));
        case 'tukey'
            if abs(r) <= delta
                sqrt_rho_prime = sqrt((1 - (r/delta)^2)^2);
            else
                sqrt_rho_prime = 0;
            end
        case 'none'
            sqrt_rho_prime = 1 ;
        otherwise
            error('Unknown loss type: %s', loss_type);
    end
    w = sqrt_rho_prime;  % optional, for debugging
end

