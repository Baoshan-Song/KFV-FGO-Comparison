classdef (Abstract) Estimator < handle
    properties
        config
        data     % Save all the data from font-end
        results
    end


    methods
        function obj = Estimator(config, data)
            obj.config = config;
            obj.data = data;
            obj.results = struct(); % All results, including intermediate and final ones
        end
    end

    methods (Abstract) 
        run(obj);
    end
    
end

