clc;clearvars;
ellipsoid = ReferenceEllipsoid('hayford');

 
% geocentric
% B = 41.27737178;
% L = 39.21259856;

% reduced 
% beta = 39.08885772;
% L = 32;

% [x, y, z] = geog2geoc(ellipsoid, beta, L, 'reduced');

%%
% 
% x = 3794779.641;
% y = 3072951.963;
% z = 4089718.725;
% 
% [b, l, h] = geoc2geog(ellipsoid, x, y, z);

%%

% JTP2:
B1 = 40.19043611;
L1 = 33.20710278;
B2 = 41.74048611;
L2 = 33.72431111;

[A1, A2, S] = JTP2(ellipsoid, B1, L1, B2, L2);