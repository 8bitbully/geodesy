%   reference:
%       Bektas, Matematik Jeodezi, p. 205
%   args:
%       ellipsoid: ReferenceEllipsoid object
%       B: geographic latitude
%       L: geographic longitude
%       L0: Longitude of the central meridian
%   returns:
%       Yg: Gauss-Kruger coordinate, Y
%       Xg: Gauss-Kruger coordinate, X
%       C: Meridian convergence
function [Yg, Xg, C] = geog2gk(ellipsoid, B, L, L0)
    l = (L - L0) * pi / 180;
    G = lat2eqdist(ellipsoid, B);
    
    eu2 = ellipsoid.SecondEccentricity^2;
    a  = ellipsoid.SemimajorAxis ;
    b  = ellipsoid.SemiminorAxis ;
    c  = a^2 / b ;
    
    sinf = @sind;
    cosf = @cosd;
    tanf = @tand;
    
    t = tanf(B);
    nu2 = eu2 * cosf(B)^2;
    V2 = 1 + nu2;
    N = c / sqrt(V2);

    C = sinf(B) * (l * 180 / pi) * ...
        (1 + 1.0153914*10e-4 * ...
        (1 + 3 * nu2) * cosf(B)^2 * (l * 180 / pi)^2);
    C = degrees2dms(C);
    
    Xg = G + ...
        (N / 24) * sinf(B) * cosf(B) * l^2 * ...
        (12 + (5 - t^2 + 9 * nu2 + 4 * nu2^2) * cosf(B)^2 * l^2);
    
    Yg = N * cosf(B) * l * ...
        (1 + 1/6 * (1 - t^2 + nu2) * cosf(B)^2 * l^2 + ...
        1/120 * (5 - 18*t^2 + t^4 - 58 * nu2 * t^2) * ...
        cosf(B)^4 * l^4);
end