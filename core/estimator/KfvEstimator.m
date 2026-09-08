classdef KfvEstimator < Estimator
    properties

    end
    
    methods
        function results = run(obj)
            % EKF derives a closed-form solution. However, iEKF, rEKF are
            % heuristic, i.e. there are not closed-form solution for them.
            % Hence, we divide them into different realization here.
            
            is_model_dataset = isobject(obj.data) && ...
                ismethod(obj.data, 'getPropagationInput') && ...
                ismethod(obj.data, 'getMeasurement');

            % Steps count and initial state
            N = obj.data.num_steps;
            if is_model_dataset
                x0 = obj.data.initial_state + obj.config.KFV.errX0;
                z = [];
                emitter = [];
            else
                z = obj.data.toa_measurements;
                emitter = obj.data.emitter_positions;
                x0 = [obj.data.true_positions(:,1); obj.data.true_velocities(:,1)] + ...
                    obj.config.KFV.errX0;
            end

            % Input parameters
            mode = obj.config.KFV.mode;
            state_size = size(x0, 1);
            if isfield(obj.config, 'state_dim') && state_size ~= obj.config.state_dim
                error('KFV:StateDimensionMismatch', ...
                    'Initial state has %d elements but config.state_dim is %d.', ...
                    state_size, obj.config.state_dim);
            end
            dt  = obj.config.KFV.dt;

            P0 = obj.config.KFV.P0;

            % System model
            f = obj.config.KFV.f;
            F = obj.config.KFV.F;
            Q = obj.config.KFV.Q; 
            omega = obj.config.KFV.omega;

            % Measuring model
            h = obj.config.KFV.h;
            H = obj.config.KFV.H;
            R = obj.config.KFV.R;
            observation_function = [];
            if isfield(obj.config.KFV, 'observation_function')
                observation_function = obj.config.KFV.observation_function;
            end

            % KF variant parameters
            max_iter = obj.config.KFV.max_iteration;
            thres_iter = obj.config.KFV.thres_iteration;
            robust_kernel = obj.config.KFV.robust_kernel;
            robust_delta = obj.config.KFV.robust_delta;
            window_size = obj.config.KFV.window_size;


            % Store all results from KFV filters
            x_est_EKF = zeros(state_size,N); x_est_iEKF = zeros(state_size,N);
            x_est_SW = zeros(state_size,N); x_est_Robust = zeros(state_size,N);
            x_est_riekf = zeros(state_size,N);

            x_est_EKF(:,1) = x0; x_est_iEKF(:,1) = x0; x_est_SW(:,1) = x0; x_est_Robust(:,1) = x0;x_est_riekf(:,1)=x0;
            P_EKF = P0; P_iEKF = P0; P_SW = repmat(P0,[1,1,N]); % window size = N
            P_Robust = P0; P_riekf = P0;

            % Only for sliding window EKF
            x_window = repmat(x0,1,window_size);
            z_window = NaN(1,window_size);

            debug_info = [];

            for k=2:N
                if is_model_dataset
                    dt_k = obj.data.getTimeStep(k);
                    process_input = obj.data.getPropagationInput(k);
                    measurement_k = obj.data.getMeasurement(k);
                    emitter_k = [];
                    R_k = [];
                else
                    dt_k = dt;
                    process_input = omega;
                    measurement_k = z(:, k);
                    emitter_k = emitter;
                    R_k = R;
                end

                % 
                % tic;
                if strcmp(mode, 'EKF')
                    % === EKF ===
                    [x_upd_EKF, P_upd_EKF, x_pred_EKF, P_pred_EKF, debug_info_EKF] = ...
                        ekf(x_est_EKF(:,k-1), P_EKF, dt_k, process_input, f, F, Q, ...
                        measurement_k, emitter_k, h, H, R_k, observation_function);
                    x_est_EKF(:,k) = x_upd_EKF; P_EKF = P_upd_EKF;
                
                % elapsed = toc;

                    debug_info{k} = debug_info_EKF;
                % debug_info{k}.time = elapsed; 
                elseif strcmp(mode, 'iEKF')
                    % === iEKF ===
                    [x_upd_iEKF, P_upd_iEKF, x_pred_iEKF, P_pred_iEKF, debug_info_iEKF] = ...
                        miekf(x_est_iEKF(:,k-1), P_iEKF, dt_k, process_input, f, F, Q, ...
                        measurement_k, emitter_k, h, H, R_k, max_iter, thres_iter, observation_function);
                    x_est_iEKF(:,k) = x_upd_iEKF; P_iEKF = P_upd_iEKF;

                    debug_info{k} = debug_info_iEKF;

                elseif strcmp(mode, 'rEKF')
                    % === Robust EKF ===
                    [x_upd_Robust, P_upd_Robust, x_pred_Robust, P_pred_Robust, debug_info_rEKF] = ...
                        rekf(x_est_Robust(:,k-1), P_Robust, dt_k, process_input, f, F, Q, ...
                        measurement_k, emitter_k, h, H, R_k, robust_kernel, robust_delta, observation_function);
                    x_est_Robust(:,k) = x_upd_Robust; P_Robust = P_upd_Robust;

                    debug_info{k} = debug_info_rEKF;

                elseif strcmp(mode, 'riEKF')
                    % === riekf ===
                    [x_upd_riEKF, P_upd_riEKF, x_pred_riEKF, P_pred_riEKF, debug_info_riEKF] = ...
                        rmiekf(x_est_riekf(:,k-1), P_riekf, dt_k, process_input, f, F, Q, ...
                        measurement_k, emitter_k, h, H, R_k, max_iter, thres_iter, ...
                        robust_kernel, robust_delta, observation_function);
                    x_est_riekf(:,k) = x_upd_riEKF; P_riekf = P_upd_riEKF;

                    debug_info{k} = debug_info_riEKF;
  
                else
                    % disp('Only support EKF/iEKF/rEKF/riEKF yet!');
                end
            end

            switch mode
                case 'EKF'
                    results.X = x_est_EKF;
                    results.debug_info = debug_info;
                case 'iEKF'
                    results.X = x_est_iEKF;
                    results.debug_info = debug_info;
                case 'rEKF'
                    results.X = x_est_Robust;
                    results.debug_info = debug_info;
                case 'riEKF'
                    results.X = x_est_riekf;
                    results.debug_info = debug_info;
                otherwise
                    disp('Only support EKF/iEKF/rEKF/riEKF yet!');
            end

        end

        function fgo_config = convert_KFV_config_to_FGO(obj)

            fgo_config.FGO.dt = obj.config.KFV.dt;
            fgo_config.FGO.errX0 = obj.config.KFV.errX0;
            fgo_config.FGO.P0 =  obj.config.KFV.P0;

            % predict
            fgo_config.FGO.f =  obj.config.KFV.f;
            fgo_config.FGO.F = obj.config.KFV.F; 
            fgo_config.FGO.Q = obj.config.KFV.Q; 
            fgo_config.FGO.omega = obj.config.KFV.omega;

            % update
            fgo_config.FGO.h =  obj.config.KFV.h; 
            fgo_config.FGO.H = obj.config.KFV.H; 

            fgo_config.FGO.R =  obj.config.KFV.R; 
            if isfield(obj.config.KFV, 'observation_function')
                fgo_config.FGO.observation_function = obj.config.KFV.observation_function;
            end
            if isfield(obj.config, 'state_dim')
                fgo_config.state_dim = obj.config.state_dim;
            end
            % variant property
            fgo_config.FGO.max_iteration  =  obj.config.KFV.max_iteration ;
            fgo_config.FGO.thres_iteration = obj.config.KFV.thres_iteration ;
            fgo_config.FGO.robust_kernel  =  obj.config.KFV.robust_kernel ;
            fgo_config.FGO.robust_delta  =  obj.config.KFV.robust_delta ;
            fgo_config.FGO.window_size    = obj.config.KFV.window_size ;

            % whether using FGO template
            fgo_config.FGO.imitate_KFV = 1;
            
            switch obj.config.KFV.mode
                case 'EKF'
                    fgo_config.FGO.max_iteration  =  1;
                    fgo_config.FGO.robust_kernel  = 'none';
                    fgo_config.FGO.window_size    = 1;
                case 'iEKF'
                    fgo_config.FGO.robust_kernel  = 'none';
                    fgo_config.FGO.window_size    = 1;
                case 'rEKF'
                    fgo_config.FGO.max_iteration  =  1;
                    fgo_config.FGO.window_size    = 1;
                case 'riEKF'
                    fgo_config.FGO.window_size    = 1;
                otherwise
                    disp('Only support EKF/iEKF/rEKF/riEKF yet!');
            end


        end
    end
end

