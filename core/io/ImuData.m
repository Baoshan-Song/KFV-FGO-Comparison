classdef ImuData < handle
    % IMU CSV parsing and GNSS-time alignment intentionally reproduce the
    % pipeline in https://github.com/ZzPolyU/KF-GNSS-INS.
    properties
        data
        batches
    end

    methods
        function obj = ImuData(imu_file)
            values = readmatrix(imu_file);
            obj.data.time = values(:, 1) * 1e-9;
            obj.data.quaternion = values(:, 5:8); % [x y z w]
            obj.data.acceleration = values(:, 30:32);
        end

        function alignToGnssTime(obj, observation_time, leap_seconds)
            if nargin < 3
                leap_seconds = 18;
            end

            imu_time = obj.data.time;
            acceleration = obj.data.acceleration;
            quaternion = obj.data.quaternion;

            [imu_time, unique_indices] = unique(imu_time, 'stable');
            acceleration = acceleration(unique_indices, :);
            quaternion = quaternion(unique_indices, :);

            for i = 2:length(imu_time)
                if imu_time(i) <= imu_time(i - 1)
                    imu_time(i) = imu_time(i - 1) + 1e-9;
                end
            end

            gps_unix_offset = seconds(datetime(1980, 1, 6) - datetime(1970, 1, 1));
            gnss_unix = zeros(observation_time.n, 1);
            for i = 1:observation_time.n
                gps_seconds = observation_time.week(i) * 7 * 24 * 3600 + observation_time.tow(i);
                gnss_unix(i) = gps_seconds + gps_unix_offset - leap_seconds;
            end

            obj.batches = cell(observation_time.n - 1, 1);
            for i = 1:observation_time.n - 1
                start_time = gnss_unix(i);
                end_time = gnss_unix(i + 1);
                indices = find(imu_time >= start_time & imu_time <= end_time);

                acceleration_start = interp1(imu_time, acceleration, start_time, 'linear');
                acceleration_end = interp1(imu_time, acceleration, end_time, 'linear');
                quaternion_start = interpolateQuaternion(imu_time, quaternion, start_time);
                quaternion_end = interpolateQuaternion(imu_time, quaternion, end_time);

                obj.batches{i} = struct( ...
                    'time', [start_time; imu_time(indices); end_time], ...
                    'acc', [acceleration_start; acceleration(indices, :); acceleration_end], ...
                    'quat', [quaternion_start; quaternion(indices, :); quaternion_end]);
            end
        end

        function batch = getBatch(obj, transition_index)
            batch = obj.batches{transition_index};
        end
    end
end

function interpolated = interpolateQuaternion(time, quaternion, query_time)
[~, closest] = min(abs(time - query_time));
if time(closest) <= query_time
    first = closest;
    second = min(closest + 1, length(time));
else
    first = max(closest - 1, 1);
    second = closest;
end

fraction = (query_time - time(first)) / (time(second) - time(first));
q1 = quaternion(first, :);
q2 = quaternion(second, :);
dot_product = dot(q1, q2);
if dot_product < 0
    q2 = -q2;
    dot_product = -dot_product;
end

if dot_product > 0.9995
    interpolated = (1 - fraction) * q1 + fraction * q2;
    interpolated = interpolated / norm(interpolated);
    return;
end

theta_0 = acos(dot_product);
sin_theta_0 = sin(theta_0);
theta = theta_0 * fraction;
interpolated = (cos(theta) - dot_product * sin(theta) / sin_theta_0) * q1 + ...
    (sin(theta) / sin_theta_0) * q2;
end
