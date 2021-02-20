%   reference:
%       Bektas, Matematik Jeodezi, p. 217
%   args:
%       Yg: Gauss-Kruger coordinate, Y
%       Xg: Gauss-Kruger coordinate, X
%       degrees: Slice degrees
%   returns:
%       right: Northing coordinates,
%       upper: Easting coordinates
function [right, upper] = gk2utm(Yg, Xg, degrees)
    % m = m0, scale factor,
    if '3' == degrees
        m = 1.0000;
    elseif '6' == degrees
        m = 0.9996;
    else
        error('Slice width incompatible.')
    end
    
    right = m * Yg + 500000;
    upper = m * Xg; % Northern hemisphere,
    % upper = m * Xg + 10000000; Southern hemisphere,
end