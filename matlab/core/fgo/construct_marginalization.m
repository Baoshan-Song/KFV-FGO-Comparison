function [estimator] = construct_marginalization(estimator, states_to_remove_gids)

%% Variables and parameters initialization
J1 = [];
J2 = [];
r1 = [];
r2 = [];

active_factors = [];
active_states = [];

% Find the factors to marginalize
for i = 1:length(estimator.factors)
    if ~strcmp(estimator.factors{i}.status, "Margin")
        estimator.factors{i} = estimator.factors{i}.evaluate();
        active_factors{end+1} = estimator.factors{i};
    end
end

% Find the states to marginalize
for i = 1:length(estimator.states)
    if ~strcmp(estimator.states(i).status, "Margin")
        active_states = [active_states,estimator.states(i)];
    end
end

% Get the number of state variables (n_x) and residuals (n_b)
n_x = size(estimator.states(1).value, 1); 
n_b = size(estimator.factors{1}.b, 1);  
m_x = 0; % Non-marginalized states
m_b = size(active_factors,2); % Number of residual factors

% Recompute non-marginalized states
for i = 1:length(estimator.states)
    if ~strcmp(estimator.states(i).status, "Margin")
        m_x = m_x + 1;
    end
end



%% Get the Jacobian and residual from latest normal equation
% (TODO) Shrink the J and r for less space and higher efficiency
J = zeros(m_b*n_b, m_x*n_x);
r = zeros(m_b*n_b,1);
row = 1; col = 1;
for i = 1:length(active_factors)
    factor = active_factors{i};
    
    % residual matrix
    sz_b = size(factor.b);
    % r((i-1)*sz_b(1)+1 : i*sz_b(1), 1 : sz_b(2)) = factor.b;
 
    r(row : row + sz_b(1)-1, 1 : sz_b(2)) = factor.b;
   
    for j = 1:length(factor.states)
        lid = factor.states(j).lid;
        % Jacobian matrix
        % J((i-1) * n_b + 1 : i*n_b, (lid-1) * n_x + 1 : lid*n_x) = ...
        %         factor.A(:, (j-1) * n_x + 1 : j*n_x);
        J(row : row + sz_b(1)-1, (lid-1) * n_x + 1 : lid*n_x) = ...
                factor.A(:, (j-1) * n_x + 1 : j*n_x);

    end
 row = row + sz_b(1);       


end


remain_states = [];
% (TODO) reorder  J/r for J1/r1/J2/r2, J1 J2 may be not at the same size
for i = 1:length(active_states)
    cur_lid = active_states(i).lid;
    if ismember(active_states(i).gid, states_to_remove_gids)

        J1 = [J1, J(:,  (cur_lid-1) * n_x + 1 : cur_lid*n_x)];
        % r1 = [r1; r((cur_lid-1)*sz_b(1)+1 : cur_lid*sz_b(1),:)];
    else
       
        J2 = [J2, J(:,  (cur_lid-1) * n_x + 1 : cur_lid*n_x)];
        remain_states = [remain_states,active_states(i)];
    end
end

tic;

%% Compute the Schur complement
H11 = J1' * J1;  % H11
H12 = J1' * J2;  % H12
H21 = J2' * J1;  % H21
H22 = J2' * J2;  % H22

b1 = J1' * r;  % b1
b2 = J2' * r;  % b2

% Check if H11 is positive definite (invertible)
try
    chol(H11);
catch
    error('H11 is singular, cannot marginalize!');
end

% Compute marginalized matrix/vector
H_marg = H22 - H21 / H11 * H12;
b_marg = b2 - H21 / H11 * b1;

% (wrong and useless) SVD decomposition for J0 and r0
% [~, S, V] = svd(H_marg);
% J0 = S.^(1/2)*V';
% r0 = S.^(1/2)*V'*b_marg;

[U, S, ~] = svd(H_marg);      % H = U*S*U'
J0 = sqrt(S) * U';        % J^T J = U*S*U' = H
r0 = (J0') \ b_marg;

% (debug) neglect the effect of marginalization factor --sbs
% if(states_to_remove_gids~=1)
% J0 = zeros(size(J0));
% r0 = zeros(size(r0));
% end


elapsed = toc;
estimator.margin_time = elapsed; 

%% Delete old state and related factors
% Delete related state
for i = 1:length(estimator.states)
    if ismember(estimator.states(i).gid, states_to_remove_gids)
        estimator.states(i).status = "Margin"; 
        estimator.states(i).lid = 0;            % Set local_id invalid 
    end
end

% Delete related factor
for i = 1:length(estimator.factors)
    for j = 1:length(estimator.factors{i}.states)
        if ismember(estimator.factors{i}.states(j).gid, states_to_remove_gids)
            estimator.factors{i}.status = "Margin"; 
            break; 
        end
    end
end

% Reassign local IDs for remaining states
lid_counter = 1; % New local_id
for i = 1:length(estimator.states)
    % if strcmp(estimator.states(i).status, "forward") % (TODO) Test for RTSS
    if  ~ismember(estimator.states(i).gid, states_to_remove_gids)
        if  ~strcmp(estimator.states(i).status, "Margin")
            estimator.states(i).lid = lid_counter;
            lid_counter = lid_counter + 1;
        end
    end
end

% add prior factor from marginalization
new_factor = margin_factor(remain_states, J0, r0, eye(size(H_marg)));
estimator = estimator.addFactor(new_factor);



end
