classdef AutoDiffRobustFactor < factor
    properties
        error_func
        kernel_func
    end

    methods
        function obj = AutoDiffRobustFactor(states, z, error_func, kernel_func)
            obj@factor(states, z, eye(1)); % Omega 在 kernel 内部处理
            obj.error_func = error_func;
            obj.kernel_func = kernel_func;
        end

        function obj = evaluate(obj)
            x = obj.states{1}.value;

            % 整个 error + kernel 组合
            total_func = @(x) obj.kernel_func(obj.error_func(x, obj.z));

            % 自动数值 Jacobian
            obj.b = total_func(x);
            obj.A = obj.numerical_jacobian(total_func, x);

        end

        function J = numerical_jacobian(~, f, x)
            fx = f(x);
            n = length(x);
            h = 1e-6;
            J = zeros(length(fx), n);
            for i = 1:n
                x1 = x; x1(i) = x1(i) + h;
                x2 = x; x2(i) = x2(i) - h;
                J(:, i) = (f(x1) - f(x2)) / (2*h);
            end
        end
    end
end
