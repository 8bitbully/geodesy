%   reference:
%       Bektas, Matematik Jeodezi, p. 217
%   args:
%       right: Northing coordinates,
%       upper: Easting coordinates
%       degrees: Slice degrees
%   returns:
%       Yg: Gauss-Kruger coordinate, Y
%       Xg: Gauss-Kruger coordinate, X
function [Yg, Xg] = utm2gk(right, upper, degrees)
    % m = m0, scale factor,
    if '3' == degrees
        m = 1.0000;
    elseif '6' == degrees
        m = 0.9996;
    else
        error('Slice width incompatible.')
    end
    
    Yg = (right - 500000) / m;
    Xg = upper / m; % Northern hemisphere,
    % Xg = (right - 10000000) / m; Southern hemisphere,
end