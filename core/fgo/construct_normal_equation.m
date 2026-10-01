function [estimator] = construct_normal_equation(estimator)

active_factors = [];

%% Relinearize the factors with given initial states
for i = 1:length(estimator.factors)
    if ~strcmp(estimator.factors{i}.status, "Margin")
        % (TODO): evaluate with latest state values --sbs
        estimator.factors{i} = estimator.factors{i}.evaluate();
        active_factors{end+1} = estimator.factors{i};
    end
end

m_x = 0;
for i = 1:length(estimator.states)
    if ~strcmp(estimator.states(i).status, "Margin")
        m_x = m_x + 1;
    end
end


%% Construct total J and r
n_x = size(estimator.states(1).value, 1); % row size of state vector
row_count = 0;
for i = 1:length(active_factors)
    row_count = row_count + size(active_factors{i}.b, 1);
end
J = zeros(row_count, m_x*n_x);
r = zeros(row_count, 1);

cur_row = 0;
for i = 1:length(active_factors)
    factor = active_factors{i};

    % residual matrix
    sz_b = size(factor.b);
    r(cur_row+1 : cur_row+sz_b(1), 1 : sz_b(2)) = factor.b;

    for j = 1:length(factor.states)
        lid = factor.states(j).lid;
        % Jacobian matrix
        J(cur_row+1 : cur_row+sz_b(1) , (lid-1) * n_x + 1 : lid*n_x) = ...
                factor.A(:, (j-1) * n_x + 1 : j*n_x);

    end
    cur_row = cur_row +sz_b(1);
end

estimator.J = J;
estimator.r = r;

end
