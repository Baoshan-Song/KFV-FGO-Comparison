classdef FgoEstimator < Estimator
    properties

    end

    methods
        function results = run(obj)

            N = obj.data.num_steps;
            range_meas = obj.data.toa_measurements;
            emitter = obj.data.emitter_positions;

            % Input parameters
            % Initial state
            x0 = [obj.data.true_positions(:,1); obj.data.true_velocities(:,1)] + obj.config.FGO.errX0; % input from data
            state_size = size(x0);
            dt  = obj.config.FGO.dt;
            P0 = obj.config.FGO.P0;

            % System model
            f = obj.config.FGO.f;
            F = obj.config.FGO.F;
            Q = obj.config.FGO.Q;

            % Measuring model
            % h = obj.config.FGO.h;
            % H = obj.config.FGO.H;
            R = obj.config.FGO.R;

            % KF variant parameters
            max_iter = obj.config.FGO.max_iteration;
            thres_iter = obj.config.FGO.thres_iteration;
            robust_kernel = obj.config.FGO.robust_kernel;
            robust_delta = obj.config.FGO.robust_delta;
            window_size = obj.config.FGO.window_size;

            % === recursive factor graph optimization (Re-FGO) and sliding window FGO(SW-FGO) ===
            %% initialize FGO estimator
            estimator = factor_graph(obj.config);

            % first node
            first_position = x0;
            first_state = state(1, 1, first_position);

            % first factor
            z = first_position; % initial position guess
            Omega = inv(P0);
            first_factor = position_factor([first_state], z, Omega);
            estimator = estimator.addState(first_state);
            estimator = estimator.addFactor(first_factor);

            for i = 2: N
                %% one-step predict
                % current state
                tic;
                current_state = estimator.states(i-1);

                % new state from propagation model
                % new_state_value = f( current_state.value , dt) ; % for constant velocity
                % new_state_value = f( current_state.value , obj.config.FGO.radius, dt) ;
                % new state from SPP
                % position0 = new_state_value(1:2);
                % [position_opt, resnorm] = solve_position_LS(position0, emitter, range_meas(:, i));
                % new_state_value(1:2,:) = position_opt(1:2);
                % new_state = state(i, i, new_state_value);    % simulate batch
                % new_state = state(i, 2, new_state_value);  % simulate recursive KF

                % initialize local id as current win_size+1
                new_state_value = f( current_state.value , dt, obj.config.FGO.omega) ; % for constant angular velocity
                new_state = state(i, estimator.win_size+1, new_state_value);  % simulate recursive KF


                % pdr propagation factor example
                % z = zeros(state_size);
                % Omega = inv(Q);
                % prop_factor = pdr_factor([current_state, new_state], z, Omega);

                % propagation factor using given functions
                prop_factor = PropagateFactor([current_state, new_state], obj.config);
                estimator = estimator.addState(new_state);
                estimator = estimator.addFactor(prop_factor);

                elapsed = toc;
                debug_info{i}.add_state_time = elapsed;

                % Estimate here for KF predict simulation
                % estimator = estimator.estimate();

                % set the first state as anchor in the SW-FGO mode
                if window_size > 1
                    if size(estimator.states,2) == 2
                        estimator = estimator.marginalize(i-1);
                        debug_info{i}.margin_time = estimator.margin_time;
                    end
                end

                % Marginalization: on -- for Re-FGO and SW-FGO, off -- for batch
                if estimator.win_size > window_size
                    % store states before marginalization for debug
                    state_before_marg=[];
                    for i=1:length(estimator.states)
                        state_before_marg=[state_before_marg,estimator.states(i).value];
                    end
                    % state_before_marg

                    % marginalize the oldest state out of window
                    estimator = estimator.marginalize(i-window_size);
                    debug_info{i}.margin_time = estimator.margin_time;
                end


                %% measurement update
                % measurement factor
                for j = 1: size(emitter, 2)
                    raw_meas.range = range_meas(j, i);
                    raw_meas.emitter = emitter(:, j);
                    raw_meas.loss_type = obj.config.FGO.robust_kernel;
                    raw_meas.loss_delta = obj.config.FGO.robust_delta;
                    raw_meas.autoDiff = obj.config.FGO.autoDiff;

                    Omega = inv(R);
                    range_factor = RangeFactor([estimator.states(i)], raw_meas, Omega);
                    estimator = estimator.addFactor(range_factor);
                end

                % estimate with Gauss-Newton method
                estimator = estimator.estimate();
                debug_info{i}.ls_time = estimator.ls_time;

                % Only for Re-FGO
                if window_size == 1
                    estimator = estimator.mar_measurements(i);
                end
                debug_info{i}.margin_meas_time = estimator.margin_meas_time;
                debug_info{i}.residual_norm_all = estimator.residual_norm_all;
            end

            % batch FGO
            % estimator = estimator.estimate();

            % output estimation results
            est_positions=[];
            for i=1:length(estimator.states)
                est_positions=[est_positions,estimator.states(i).value];
            end

            results.X = est_positions;
            results.debug_info = debug_info;

            % single point positioning (SPP)
            function [position_opt, resnorm] = solve_position_LS(position0, emitter, range_meas, max_iter)
                if nargin < 4
                    max_iter = 10;
                end

                position = position0(:);  % Ensure column vector

                for iter = 1:max_iter
                    r = zeros(size(range_meas));
                    J = zeros(length(range_meas), 2);

                    for j = 1:length(range_meas)
                        e_j = emitter(:, j);
                        diff = position - e_j;
                        dist = norm(diff);
                        r(j) = dist - range_meas(j);

                        if dist > 1e-6
                            J(j, :) = (diff / dist)';
                        else
                            J(j, :) = [0, 0]; % Avoid NaN
                        end
                    end

                    % Gauss-Newton step
                    delta = -(J' * J) \ (J' * r);
                    position = position + delta;

                    if norm(delta) < 1e-6
                        break;
                    end
                end

                position_opt = position;
                resnorm = norm(r)^2;
            end

        end


    end

end

