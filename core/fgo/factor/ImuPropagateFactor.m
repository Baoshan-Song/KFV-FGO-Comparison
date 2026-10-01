classdef ImuPropagateFactor < factor
    properties
        propagation_input
        dt
    end

    methods
        function obj = ImuPropagateFactor(states, propagation_input, dt, config)
            obj@factor(states, config);
            obj.propagation_input = propagation_input;
            obj.dt = dt;
        end

        function obj = evaluate(obj)
            current_state = obj.states(1);
            next_state = obj.states(2);
            predicted = obj.config.FGO.f(current_state.value, obj.dt, obj.propagation_input);
            F = obj.config.FGO.F(current_state.value, obj.dt, obj.propagation_input);
            state_size = length(next_state.value);

            obj.A = [F, -eye(state_size)];
            obj.b = next_state.value - predicted;
            obj.Omega = inv(obj.config.FGO.Q);
            square_root_information = chol(obj.Omega, 'lower');
            obj.A = square_root_information' * obj.A;
            obj.b = square_root_information' * obj.b;
        end
    end
end
