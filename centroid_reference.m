function cref = centroid_reference(type, varargin)
%CENTROID_REFERENCE Factory for centroid reference trajectories
%
% Usage:
%   cref = centroid_reference("circle")
%   cref = centroid_reference("ellipse")
%   cref = centroid_reference("static",[5;5])
%
% Output:
%   cref(t) -> 2x1 centroid position

    switch lower(type)

        case "ellipse"
            % default elliptical reference
            omega = 0.1;
            ax = 2;
            ay = 3;
            offset = 7;

            cref = @(t) [ ...
                ax*sin(omega*t); ...
                ay*cos(omega*t) ] + offset;

        case "circle"
            omega = 0.1;
            r = 3;
            offset = [7;7];

            cref = @(t) offset + r * [ ...
                cos(omega*t); ...
                sin(omega*t) ];

        case "static"
            p = varargin{1};   % 2x1 vector
            cref = @(t) p;

        case "line"
            if isempty(varargin)
                p0 = [3; 0];
                v  = [0.1; -0.4];
            else
                p0 = varargin{1};
                v  = varargin{2};
            end
        
            cref = @(t) p0 + v*t;

        case "lemniscate"
            omega = 0.1;
            a = 3;
            offset = 7;

            cref = @(t) offset + (a*sqrt(2)) * [ ...
                cos(omega*t) ./ (1 + sin(omega*t).^2); ...
                sin(omega*t).*cos(omega*t) ./ (1 + sin(omega*t).^2) ];

        otherwise
            error("Unknown centroid reference type: %s", type)
    end
end
