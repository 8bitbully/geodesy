%   args:
%       ellipsoid: ReferenceEllipsoid object
%       B: geographic latitude
%       L: geographic longitude
%   returns:
%       Yg: Gauss-Kruger coordinate, Y
%       Xg: Gauss-Kruger coordinate, X
%       C: Meridian convergence
function [Yg, Xg, C] = geographic2GaussKruger(ellipsoid, B, L)
    L0 = src.dom(L, '3');
    [Yg, Xg, C] = src.geog2gk(ellipsoid, B, L, L0);
end