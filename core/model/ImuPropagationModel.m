classdef ImuPropagationModel < handle
    % Simplified propagation intentionally reproduces KF-GNSS-INS: use the
    % first valid sample and external attitude, then apply first-order Euler
    % position/velocity integration without an additional 0.5*a*dt^2 term.
    properties
        Q
        gravity_norm
    end

    methods
        function obj = ImuPropagationModel(Q, gravity_norm)
            obj.Q = Q;
            obj.gravity_norm = gravity_norm;
        end

        function x_pred = stateTransition(obj, x_prev, ~, imu_measurement)
            [x_pred, ~] = obj.computePropagation(x_prev, imu_measurement);
        end

        function F = transitionJacobian(obj, x_prev, ~, imu_measurement)
            [~, F] = obj.computePropagation(x_prev, imu_measurement);
        end

        function [x_pred, F, Q] = computePropagation(obj, x_prev, imu_measurement)
            if length(x_prev) ~= 10
                error('KFV:InvalidGnssInsState', 'GNSS/IMU propagation requires a 10-state vector.');
            end

            position = x_prev(1:3);
            velocity = x_prev(4:6);
            accelerometer_bias = x_prev(7:9);
            receiver_clock = x_prev(10);
            dt = imu_measurement.time(end) - imu_measurement.time(1);

            sample = 1;
            while sample <= size(imu_measurement.acc, 1) && ...
                    (any(isnan(imu_measurement.acc(sample, :))) || ...
                    any(isnan(imu_measurement.quat(sample, :))))
                sample = sample + 1;
            end
            if sample > size(imu_measurement.acc, 1)
                error('KFV:InvalidImuBatch', 'IMU batch contains no valid acceleration/quaternion pair.');
            end

            acceleration_body = imu_measurement.acc(sample, :)' - accelerometer_bias;
            quaternion = imu_measurement.quat(sample, :); % [x y z w]
            rotation_navigation_body = quaternionToRotationMatrix( ...
                [quaternion(4), quaternion(1), quaternion(2), quaternion(3)]);
            rotation_ecef_navigation = obj.localToEcefRotation(position);
            rotation_ecef_body = rotation_ecef_navigation * rotation_navigation_body;
            gravity_navigation = [0; 0; -obj.gravity_norm];
            acceleration_ecef = rotation_ecef_body * acceleration_body + ...
                rotation_ecef_navigation * gravity_navigation;

            % Intentionally matches the source repository's propagation.
            position = position + velocity * dt;
            velocity = velocity + acceleration_ecef * dt;
            x_pred = [position; velocity; accelerometer_bias; receiver_clock];

            F = eye(10);
            F(1:3, 4:6) = eye(3) * dt;
            F(4:6, 7:9) = -rotation_ecef_body * dt;
            Q = obj.Q;
        end
    end

    methods (Static, Access = private)
        function rotation = localToEcefRotation(position_ecef)
            position = gt.Gpos(position_ecef', 'xyz');
            longitude = deg2rad(position.llh(2));
            latitude = deg2rad(position.llh(1));
            rotation = [ ...
                -sin(longitude), -sin(latitude) * cos(longitude),  cos(latitude) * cos(longitude); ...
                 cos(longitude), -sin(latitude) * sin(longitude),  cos(latitude) * sin(longitude); ...
                 0,               cos(latitude),                   sin(latitude)];
        end
    end
end

function rotation = quaternionToRotationMatrix(quaternion)
quaternion = quaternion(:)' / norm(quaternion);
w = quaternion(1); x = quaternion(2); y = quaternion(3); z = quaternion(4);
rotation = [ ...
    1 - 2 * (y^2 + z^2), 2 * (x*y - z*w),     2 * (x*z + y*w); ...
    2 * (x*y + z*w),     1 - 2 * (x^2 + z^2), 2 * (y*z - x*w); ...
    2 * (x*z - y*w),     2 * (y*z + x*w),     1 - 2 * (x^2 + y^2)];
end
