classdef PropagateFactor < factor
    properties

    end
    
    methods
        function obj = PropagateFactor(states, config)
            obj@factor(states, config);
        end

        function obj = evaluate(obj)
            % Extract states
            cur_state = obj.states(1);
            new_state = obj.states(2);
            % F = obj.config.FGO.F(cur_state.value, obj.config.FGO.dt);
             % test constant circle velocity --sbs
            % F = obj.config.FGO.F(cur_state.value, obj.config.FGO.radius, obj.config.FGO.dt);
            % test constant angular velocity --sbs
            F = obj.config.FGO.F(cur_state.value, obj.config.FGO.dt, obj.config.FGO.omega);

            % Propagation factor
            obj.A = [ F, -eye(4)];
            
            % predictde_new_state_value = obj.config.FGO.f(cur_state.value, ...
            %                                             obj.config.FGO.dt);       
             % test constant circle velocity --sbs
            % predictde_new_state_value = obj.config.FGO.f(cur_state.value, ...
            %                                             obj.config.FGO.radius, ...
            %                                             obj.config.FGO.dt);       

             predictde_new_state_value = obj.config.FGO.f( cur_state.value , ...
                 obj.config.FGO.dt, obj.config.FGO.omega) ; % test constant angular velocity --sbs

            obj.b = - (new_state.value - predictde_new_state_value) ;

            obj.Omega = inv(obj.config.FGO.Q);
            
            % (TODO): Mutiply the obj.Omega to A
            L = chol(obj.Omega, 'lower');
            obj.A = L' * obj.A;            
            obj.b = L' * obj.b;

        end
    end
end
