classdef margin_factor < factor
    properties

    end
    
    methods
        function obj = margin_factor(states, A, b, Omega)
            obj@factor(states,A,b,Omega);
        end

        function obj = evaluate(obj)
            % (TODO) Move the construct_marginalization to here
        end
    end
end
