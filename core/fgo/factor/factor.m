classdef (Abstract) factor < handle
    properties
        states    % States (only global id of states are used in factors)
        z      % Raw measurement        
        A        % Jacobian
        b        % residual
        Omega    % Information matrix

        status  % 'Margined' or ''

        auto_diff % Whether using auto numerical Jacobian by gradient(y,x)
    
        config
    end
    
    methods
        function obj = factor(varargin)  
            if nargin == 2   %(states, config)
                obj.states = varargin{1};
                obj.config = varargin{2};
            elseif nargin == 3   %(states, z, Omega)
                obj.states = varargin{1};
                obj.z = varargin{2};
                obj.Omega = varargin{3};
            elseif nargin == 4   %(states, A, b, Omega)
                obj.states = varargin{1};
                obj.A= varargin{2};
                obj.b = varargin{3};
                obj.Omega = varargin{4};
            else
                err('Too many input arguments');
            end
        end

        function obj = evaluate(obj)
            disp('This is a virtual evaluate function.');
        end
    end
end
