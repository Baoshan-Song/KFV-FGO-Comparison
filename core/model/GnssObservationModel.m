classdef GnssObservationModel < handle
    properties
        minimum_snr
        minimum_elevation_deg
        Tref
        a
        A
        Fref
    end

    methods
        function obj = GnssObservationModel(config)
            obj.minimum_snr = config.minimum_snr;
            obj.minimum_elevation_deg = config.minimum_elevation_deg;
            obj.Tref = config.Tref;
            obj.a = config.a;
            obj.A = config.A;
            obj.Fref = config.Fref;
        end

        function [residual, H, R] = computeObservation(obj, state, measurement)
            if length(state) ~= 10
                error('KFV:InvalidGnssInsState', 'GNSS observation requires a 10-state vector.');
            end

            % The same state-independent measurement packet is reused by
            % every EKF/FGO iteration. All receiver-dependent RTKLIB work
            % starts from the current state passed to this function.
            % Gsat is a handle class. Copy it so relinearization has no
            % hidden mutable state and separate estimators cannot interfere.
            satellite = measurement.satellite.copy();
            satellite.setRcvPos(gt.Gpos(state(1:3)', 'xyz'));
            observation = measurement.observation.residuals(satellite);
            pseudorange_residual = observation.L1.resPc - state(10) - measurement.tgd;

            valid = ~isnan(pseudorange_residual);
            if ~isempty(obj.minimum_elevation_deg)
                valid = valid & satellite.el >= obj.minimum_elevation_deg;
            end
            if ~isempty(obj.minimum_snr)
                valid = valid & observation.L1.S >= obj.minimum_snr;
            end

            measurement_count = nnz(valid);
            if measurement_count < 4
                error('KFV:InsufficientSatellites', ...
                    'GNSS update requires at least four valid satellites; found %d.', measurement_count);
            end

            residual = pseudorange_residual(valid)';
            H = zeros(measurement_count, length(state));
            H(:, 1) = -satellite.ex(valid)';
            H(:, 2) = -satellite.ey(valid)';
            H(:, 3) = -satellite.ez(valid)';
            H(:, 10) = 1;
            R = obj.computeCovariance(satellite.el(valid), observation.L1.S(valid));
        end
    end

    methods (Access = private)
        function R = computeCovariance(obj, elevation, snr)
            measurement_count = length(elevation);
            variances = zeros(measurement_count, 1);
            for i = 1:measurement_count
                elevation_rad = deg2rad(elevation(i));
                cn0 = snr(i);
                if ~isfinite(elevation_rad), elevation_rad = 1e-3; end
                if ~isfinite(cn0), cn0 = obj.Tref; end

                elevation_term = 1 / sin(max(elevation_rad, 1e-3))^2;
                snr_term = 10^(-(cn0 - obj.Tref) / obj.a);
                kappa = obj.A * 10^((obj.Fref - obj.Tref) / obj.a) - 1;
                correction_term = kappa * (cn0 - obj.Tref) / ...
                    (obj.Fref - obj.Tref) + 1;
                variances(i) = elevation_term * snr_term * correction_term;
            end
            R = diag(variances);
        end
    end
end
