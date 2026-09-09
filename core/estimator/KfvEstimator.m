classdef KfvEstimator < Estimator
    properties

    end
    
    methods
        function results = run(obj)
            % EKF derives a closed-form solution. However, iEKF, rEKF are
            % heuristic, i.e. there are not closed-form solution for them.
            % Hence, we divide them into different realization here.
            
            % Steps count 
            N = obj.data.num_steps;
            z = obj.data.toa_measurements;
            emitter = obj.data.emitter_positions;

            % Input parameters
            mode = obj.config.KFV.mode;
            % Initial state
            % x0 = obj.config.KFV.X0; % input from config
            x0 = [obj.data.true_positions(:,1); obj.data.true_velocities(:,1)]+ obj.config.KFV.errX0; % test different initial point --sbs
            state_size = size(x0, 1);
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
                % 
                % tic;
                if strcmp(mode, 'EKF')
                    % === EKF ===
                    [x_upd_EKF, P_upd_EKF, x_pred_EKF, P_pred_EKF, debug_info_EKF] = ...
                        ekf(x_est_EKF(:,k-1), P_EKF, dt, omega, f, F, Q, z(:, k), emitter, h, H, R);
                    x_est_EKF(:,k) = x_upd_EKF; P_EKF = P_upd_EKF;
                
                % elapsed = toc;

                    debug_info{k} = debug_info_EKF;
                % debug_info{k}.time = elapsed; 
                elseif strcmp(mode, 'iEKF')
                    % === iEKF ===
                    [x_upd_iEKF, P_upd_iEKF, x_pred_iEKF, P_pred_iEKF, debug_info_iEKF] = ...
                        miekf(x_est_iEKF(:,k-1), P_iEKF, dt, omega, f, F, Q, z(:, k), emitter, h, H, R, max_iter,thres_iter);
                    x_est_iEKF(:,k) = x_upd_iEKF; P_iEKF = P_upd_iEKF;

                    debug_info{k} = debug_info_iEKF;

                elseif strcmp(mode, 'rEKF')
                    % === Robust EKF ===
                    [x_upd_Robust, P_upd_Robust, x_pred_Robust, P_pred_Robust, debug_info_rEKF] = ...
                        rekf(x_est_Robust(:,k-1), P_Robust, dt, omega, f, F, Q, z(:, k), emitter, h, H, R, robust_kernel, robust_delta);
                    x_est_Robust(:,k) = x_upd_Robust; P_Robust = P_upd_Robust;

                    debug_info{k} = debug_info_rEKF;

                elseif strcmp(mode, 'riEKF')
                    % === riekf ===
                    [x_upd_riEKF, P_upd_riEKF, x_pred_riEKF, P_pred_riEKF, debug_info_riEKF] = ...
                        rmiekf(x_est_riekf(:,k-1), P_riekf, dt, omega, f, F, Q, z(:, k), emitter, h, H, R, max_iter,thres_iter,robust_kernel, robust_delta);
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

