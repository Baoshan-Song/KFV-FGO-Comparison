function data = validate_comparison_data(data)
%VALIDATE_COMPARISON_DATA Validate data required by KFV and FGO estimators.

if ~isstruct(data) || ~isscalar(data)
    invalid('Dataset must be a scalar struct loaded from a MAT file.');
end

required = {'true_positions', 'true_velocities', 'toa_measurements', ...
    'emitter_positions', 'num_steps'};
missing = required(~isfield(data, required));
if ~isempty(missing)
    invalid('Dataset is missing required fields: %s.', strjoin(missing, ', '));
end

n = data.num_steps;
if ~isnumeric(n) || ~isscalar(n) || ~isfinite(n) || n < 2 || fix(n) ~= n
    invalid('num_steps must be a finite integer greater than or equal to 2.');
end

validateMatrix(data.true_positions, 'true_positions', 2, n);
validateMatrix(data.true_velocities, 'true_velocities', 2, n);

toa = data.toa_measurements;
if ~isnumeric(toa) || ~ismatrix(toa) || isempty(toa) || ...
        size(toa, 2) ~= n || any(~isfinite(toa), 'all')
    invalid('toa_measurements must be a finite M-by-num_steps numeric matrix.');
end

emitters = data.emitter_positions;
if ~isnumeric(emitters) || ~ismatrix(emitters) || size(emitters, 1) ~= 2 || ...
        size(emitters, 2) ~= size(toa, 1) || any(~isfinite(emitters), 'all')
    invalid(['emitter_positions must be a finite 2-by-M numeric matrix, ' ...
        'where M matches the rows of toa_measurements.']);
end
end

function validateMatrix(value, name, rows, columns)
if ~isnumeric(value) || ~ismatrix(value) || ...
        ~isequal(size(value), [rows, columns]) || any(~isfinite(value), 'all')
    invalid('%s must be a finite %d-by-num_steps numeric matrix.', name, rows);
end
end

function invalid(message, varargin)
error('KFVFGO:InvalidData', message, varargin{:});
end
