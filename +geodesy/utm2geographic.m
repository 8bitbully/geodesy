%   reference:
%       Bektas, Matematik Jeodezi, p. 217
%   args:
%       ellipsoid: ReferenceEllipsoid object
%       right: Northing coordinates,
%       upper: Easting coordinates
%       L0: Longitude of the central meridian
%       degrees: Slice degrees
%   returns:
%       B: geographic latitude
%       L: geographic longitude
function [B, L] = utm2geographic(ellipsoid, right, upper, L0, degrees)
    [Yg, Xg] = src.utm2gk(upper, right, degrees);
    [B, L, ~] = src.gk2geog(ellipsoid, Yg, Xg, L0);
end