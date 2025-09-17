classdef factor_graph
    properties
        states     % Varaible states
        factors   % Measurement factors
        marginal_factor % Marginalization factor

        mode     % Sliding window, batch

        J             % Full Jacobian matrix
        r             % Full residual matrix
        win_size % sliding window size

        config

        latest_information_matrix % Store information matrix to marginalize used measurements (used in Re-FGO mode)
        ls_time
        margin_time
        margin_meas_time

        residual_norm_all
    end
    
    methods
        function obj = factor_graph(config)
            obj.states = [];
            obj.factors = [];

            obj.config = config;
            obj.win_size = 0;
        end
        
        function obj = addState(obj, node)
            obj.states = [obj.states, node];
            obj.win_size  = obj.win_size +1;
        end
        
        function obj = addFactor(obj, factor)
            obj.factors{end+1} = factor;
        end
        
        function obj = normal_equation(obj)
            obj = construct_normal_equation(obj);
        end

        function obj = marginalize(obj, node_to_remove)
            obj = construct_marginalization(obj, node_to_remove);
            obj.win_size = obj.win_size - size(node_to_remove,1) ;
        end

        function obj = mar_measurements(obj, related_states)
            obj = construct_meas_prior_factor(obj,related_states);
        end
        
        function obj = estimate(obj)
            obj = estimate_fgo(obj);
        end
    end
end
