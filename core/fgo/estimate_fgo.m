function estimator = estimate_fgo(estimator)

% Initialize parameters
max_iter = estimator.config.FGO.max_iteration;
tolerance = estimator.config.FGO.thres_iteration;
% lambda = 1e-3;  % Levenberg-Marquardt parameter

% Row size of state vector
n_x = size(estimator.states(1).value, 1);
m_x = 0;
states_gid = [];
for i = 1:length(estimator.states)
    if ~strcmp(estimator.states(i).status, "Margin")
        m_x = m_x + 1;
        states_gid = [states_gid; estimator.states(i).gid];
    end
end


estimator.residual_norm_all=[];
estimator.ls_time =0;

% Iteration for optimization
for iter = 1:max_iter
    tic;

    %% construct normal equation with latest initial guess
    estimator = estimator.normal_equation();
    J = estimator.J;
    r = estimator.r;

    % test: only get current state for comparison with KFV--sbs
    estimator.residual_norm_all = [estimator.residual_norm_all, [estimator.states(states_gid(1)).value; norm(r)] ];

    H = J'*J;
    % Estimate with Gaussian-Newton based methods
    x_opt = (J' * J)\ (J' * r);

    elapsed = toc;
    estimator.ls_time = estimator.ls_time +elapsed;

    % Update active states in the window
    for i = 1:m_x
        % disp('before: ');estimator.states(states_gid(i)).value
        % disp('delta_x: '); x_opt((i-1)*n_x+1 : i*n_x , 1)
        estimator.states(states_gid(i)).value = estimator.states(states_gid(i)).value + x_opt((i-1)*n_x+1 : i*n_x , 1);
        % disp('after: '); estimator.states(states_gid(i)).value
    end


    % Store the information matrix of current estimation
    estimator.latest_information_matrix = J'*J;

    if norm(x_opt)/length(x_opt)<tolerance
        break;
    end

end


end

