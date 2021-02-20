%   args:
%       ellipsoid: ReferenceEllipsoid object
%       Yg: Gauss-Kruger coordinate, Y
%       Xg: Gauss-Kruger coordinate, X
%       degrees: Slice degrees
%   returns:
%       right: Northing coordinates,
%       upper: Easting coordinates
function [right, upper] = geographic2UTM(ellipsoid, B, L, degrees) % L0
    % L0 = dom;
    L0 = src.dom(L, degrees);
    [Yg, Xg, ~] = src.geog2gk(ellipsoid, B, L, L0);
    [right, upper] = src.gk2utm(Yg, Xg, degrees);
end