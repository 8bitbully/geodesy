%   reference:
%       Bektas, Matematik Jeodezi, p. 210
%   args:
%       ellipsoid: ReferenceEllipsoid object
%       Yg: Gauss-Kruger coordinate, Y
%       Xg: Gauss-Kruger coordinate, X
%       L0: Longitude of the central meridian
%   returns:
%       B: geographic latitude
%       L: geographic longitude
%       C: Meridian convergence
function [B, L, C] = gk2geog(ellipsoid, Yg, Xg, L0)
    eu2 = ellipsoid.SecondEccentricity^2;
    a  = ellipsoid.SemimajorAxis ;
    b  = ellipsoid.SemiminorAxis ;
    c  = a^2 / b ;
    
    cosf = @cosd;
    tanf = @tand;
    
    Bf = dist2eqlat(ellipsoid, Xg); %G = Xg;
    
    ro = 180 / pi;
    t = tanf(Bf);
    nu2 = eu2 * cosf(Bf)^2;
    V2 = 1 + nu2;
    Nf = c / sqrt(V2);

    B = Bf - ...
        ro / 2 * tanf(Bf) * (Yg / Nf)^2 * ...
        (1 + nu2 - 1/12 * (5 + 6*nu2 + 3*t^2 - 6*nu2*t^2 * (Yg/Nf)^2));
    
    L = L0 + ...
        ro / cosf(B) * (Yg/Nf) * ...
        (1 - 1/6 *(1 + nu2 + 2*t^2) * (Yg/Nf)^2 + ...
        1/120 * (5 + 28*t^2 + 24*t^4) * (Yg/Nf)^4);

    C = ro * t / Nf * Yg - ...
        ro * t * ((1 + t^2 - nu2) / 3 ) * (Yg / Nf)^3;
end