clearvars; clc;
main = localfunctions;
e = ReferenceEllipsoid('hayford');

pf = Pafta.PAFTA_25;

[cidx, area, edge] = paftaEdge(pf);
[B, L] = paftaParser(pf, 39.5, 40.5);

% Pafta köşe koordinatları:
PAFTA_KOSE = [double(cidx)', B(cidx)', L(cidx)'];

% Kenar uzunlukları:
for i = 1 : size(edge, 1)
    G(i,1) = abs(lat2eqdist(e, B(edge(i, 1))) - lat2eqdist(e, B(edge(i, 2))));
    if B(edge(i, 1)) == B(edge(i, 2))
        Sp(i, 1) = longt2dist(e, L(edge(i, 1)), L(edge(i, 2)), B(edge(i, 1)));
        [~, ~, S(i, 1)] = JTP2(e, B(edge(i, 1)), B(edge(i, 2)), L(edge(i, 1)), L(edge(i, 2)));
        diff(i, 1) = Sp(i, 1) - S(i, 1);
    end
end

KENAR = [double(cidx)' ,G, Sp, S, diff];

% Pafta Alanı
for i = 1 : size(area,1)
    F(i, 1) = ellipsoidArea(e, B(area(i, 1)), B(area(i, 3)), L(area(i, 1)), L(area(i, 3)));
    [km(i, 1), hektar(i, 1), donum(i, 1)] = areaConvert(F(i, 1));
end

ALAN = [F, km, hektar, donum];

% Ağrılık merkezi yarıçapları
for i = 1: size(area, 1)
    Bort(i, 1) = (B(area(i, 1)) + B(area(i, 3))) / 2;
    [N(i, 1), M(i, 1), R(i, 1)] = radiusCurveture(e, Bort(i, 1));
end

AGIRLIK = [N, M, R];


% Jeosantrik koordinatlar
[x, y, z] = geog2geoc(e, B(cidx), L(cidx));
JEOCENTRIC = [double(cidx)', x', y', z'];

% diary hayford-25.txt
% PAFTA_KOSE
% fprintf('--------------------------------')
% KENAR
% fprintf('--------------------------------')
% ALAN
% fprintf('--------------------------------')
% AGIRLIK
% fprintf('--------------------------------')
% JEOCENTRIC
% diary off


function [km, hektar, donum] = areaConvert(m)
    km = m / 1e6;
    hektar = m / 1e4;
    donum =  m / 1e3;
end