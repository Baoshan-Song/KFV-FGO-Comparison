function [x, status] = single_point_position(observation, navigation)
%SINGLE_POINT_POSITION Initialize ECEF position and receiver clock bias.

x = zeros(4, 1);
status = false;
satellite = gt.Gsat(observation, navigation);

for iteration = 1:10
    satellite.setRcvPos(gt.Gpos(x(1:3)', 'xyz'));
    corrected_observation = observation.residuals(satellite);
    residual = corrected_observation.L1.resPc - x(4) - navigation.getTGD(satellite.sat);

    valid = ~isnan(residual) & satellite.el > 10;
    observation_count = nnz(valid);
    if observation_count < 4
        error('KFV:InsufficientSatellites', ...
            'SPP requires at least four valid satellites; found %d.', observation_count);
    end

    variance_at_90_degrees = 0.5^2;
    weights = 1 ./ (variance_at_90_degrees ./ sind(satellite.el(valid)));
    H = zeros(observation_count, 4);
    H(:, 1) = -satellite.ex(valid)';
    H(:, 2) = -satellite.ey(valid)';
    H(:, 3) = -satellite.ez(valid)';
    H(:, 4) = 1;

    correction = lscov(H, residual(valid)', weights);
    x = x + correction;
    if norm(correction) < 1e-3
        status = true;
        break;
    end
end

if ~status
    warning('KFV:SppDidNotConverge', 'SPP did not converge within 10 iterations.');
end
end
