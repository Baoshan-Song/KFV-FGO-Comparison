classdef GnssImuDataset < handle
    properties
        type = 'gnss_imu'
        num_steps
        initial_state
    end

    properties (Access = private)
        observations
        navigation
        satellite_cache
        tgd
        imu
    end

    methods
        function obj = GnssImuDataset(data_config)
            required = {'imu_file', 'obs_file', 'nav_file'};
            for i = 1:length(required)
                field_name = required{i};
                if ~isfield(data_config, field_name)
                    error('KFV:MissingDataPath', ...
                        'Missing data configuration field: %s', field_name);
                end
                if ~isfile(data_config.(field_name))
                    error('KFV:MissingDataFile', ...
                        'Data file does not exist: %s', data_config.(field_name));
                end
            end

            obj.observations = gt.Gobs(data_config.obs_file);
            obj.navigation = gt.Gnav(data_config.nav_file);
            % Make frequency/wavelength initialization explicit. The
            % original pipeline relied on a shared-handle side effect in
            % Gsat.setRcvPos(), which is lost when satellite data is cached.
            obj.observations.setFrequencyFromNav(obj.navigation);

            % Satellite position, velocity and clock are independent of the
            % estimated receiver state and are therefore cached once.
            obj.satellite_cache = gt.Gsat(obj.observations, obj.navigation);
            obj.tgd = obj.navigation.getTGD(obj.satellite_cache.sat);
            obj.num_steps = obj.observations.n;

            obj.imu = ImuData(data_config.imu_file);
            leap_seconds = 18;
            if isfield(data_config, 'leap_seconds')
                leap_seconds = data_config.leap_seconds;
            end
            obj.imu.alignToGnssTime(obj.observations.time, leap_seconds);

            spp_state = single_point_position(obj.observations.selectTime(1), obj.navigation);
            obj.initial_state = [spp_state(1:3); zeros(3, 1); zeros(3, 1); spp_state(4)];
        end

        function input = getPropagationInput(obj, epoch)
            input = obj.imu.getBatch(epoch - 1);
        end

        function dt = getTimeStep(obj, epoch)
            input = obj.getPropagationInput(epoch);
            dt = input.time(end) - input.time(1);
        end

        function measurement = getMeasurement(obj, epoch)
            % Return only state-independent data for one GNSS epoch. The
            % receiver state is intentionally not an input here: range,
            % LOS, elevation, atmosphere corrections, residuals, and R are
            % recomputed by GnssObservationModel at each linearization.
            if ~isscalar(epoch) || epoch ~= fix(epoch) || ...
                    epoch < 1 || epoch > obj.num_steps
                error('KFV:InvalidMeasurementEpoch', ...
                    'Measurement epoch must be an integer from 1 to %d.', obj.num_steps);
            end

            measurement.epoch = epoch;
            measurement.observation = obj.observations.selectTime(epoch);
            measurement.satellite = obj.satellite_cache.selectTime(epoch);
            measurement.tgd = obj.tgd;
        end
    end
end
