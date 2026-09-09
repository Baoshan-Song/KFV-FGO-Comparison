function [estimator] = construct_meas_prior_factor(estimator, states_to_remove_gids)

% Mark factors related to states_to_remove_gids as "Margin"
for i = 1:length(estimator.factors)
    factor = estimator.factors{i};
    
    % Check if this factor depends on any of the states to be removed
    for j = 1:length(factor.states)
        if ismember(factor.states(j).gid, states_to_remove_gids)
            estimator.factors{i}.status = "Margin";  % Mark factor for marginalization
            break;      % No need to continue checking once it's marked
        end
    end
end

tic;

% Add prior factor from marginalization
new_factor = position_factor(estimator.states(states_to_remove_gids), ...
    estimator.states(states_to_remove_gids).value, ...
    estimator.latest_information_matrix);
estimator = estimator.addFactor(new_factor);


elapsed = toc;
estimator.margin_meas_time = elapsed; 


end
