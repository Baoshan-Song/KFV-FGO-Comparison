classdef position_factor < factor
    properties

    end
    
    methods
        function obj = position_factor(states, z, Omega)
            obj@factor(states,z,Omega);
        end

        function obj = evaluate(obj)
            % Extract states
            state = obj.states(1);

            % Position factor
            obj.A = -eye(4);
            obj.b = obj.z -  state.value;

            
            % (TODO): Mutiply the obj.Omega to A
            L = chol(obj.Omega, 'lower');
            obj.A = L' * obj.A;            
            obj.b = L' * obj.b;

        end
    end
end
