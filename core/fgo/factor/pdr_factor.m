classdef pdr_factor < factor
    properties

    end
    
    methods
        function obj = pdr_factor(states, z, Omega)
            obj@factor(states,z,Omega);
        end

        function obj = evaluate(obj)
            % Extract states
            cur_state = obj.states(1);
            new_state = obj.states(2);

            % Propagation factor
            obj.A = [eye(4), -eye(4)];
            obj.b = obj.z - (new_state.value - cur_state.value) ;
            
            % (TODO): Mutiply the obj.Omega to A
            L = chol(obj.Omega, 'lower');
            obj.A = L' * obj.A;            
            obj.b = L' * obj.b;

        end
    end
end
