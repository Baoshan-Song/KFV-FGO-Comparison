classdef AutoDiffFactor < factor
    properties
        error_func % function handle @(x) → residual
        delta      % precision for perturbation
    end

    methods
        function obj = AutoDiffFactor(states, z, Omega, error_func, delta)
            obj@factor(states, z, Omega);
            obj.error_func = error_func;
            if nargin < 5, obj.delta = 1e-4; else, obj.delta = delta; end
        end

        function obj = evaluate(obj)
            % Get current state value
            x = obj.states(1).value;

            % Residual
            r = obj.error_func(x, obj.z);   % residual, m x 1
            J = obj.numerical_jacobian(x, obj.z, obj.error_func);  % Jacobian, m x n

            % Construct A and b
            L = chol(obj.Omega, 'lower');
            obj.b = L' * r;
            obj.A = L' * J;
        end

        function J = numerical_jacobian(obj, x, z, f)
            n = length(x);
            r0 = f(x, z);
            m = length(r0);
            J = zeros(m, n);
            h = obj.delta;

            for i = 1:n
                x1 = x; x1(i) = x1(i) + h;
                x2 = x; x2(i) = x2(i) - h;
                J(:,i) = (f(x1,z) - f(x2,z)) / (2*h);
            end
        end
    end
end
