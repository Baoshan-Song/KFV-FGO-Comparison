classdef MockModelDataset < handle
    properties
        type = 'mock_model'
        num_steps = 4
        initial_state = zeros(10, 1)
    end

    methods
        function input = getPropagationInput(~, epoch)
            input = epoch;
        end

        function dt = getTimeStep(~, ~)
            dt = 1;
        end

        function measurement = getMeasurement(~, epoch)
            measurement.H = zeros(4, 10);
            measurement.H(:, 1:3) = [ ...
                 1,  0,  0; ...
                 0,  1,  0; ...
                 0,  0,  1; ...
                -1, -1, -1];
            measurement.H(:, 10) = 1;
            measurement.z = [epoch; -epoch; 0.5 * epoch; 2 * epoch];
            measurement.R = diag([1, 2, 3, 4]);
        end
    end
end
