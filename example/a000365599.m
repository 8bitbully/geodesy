clearvars; clc;
main = localfunctions;
% elipsoid secimi
e = ReferenceEllipsoid('hayford');

% pafta secimi
pf = Pafta.PAFTA_25;

% pafta indexleri
[cidx, area, edge] = paftaEdge(pf);
[B, L] = paftaParser(pf, 39.5, 37);

% Pafta köse koordinatlari:
PAFTA_KOSE = [double(cidx)', B(cidx)', L(cidx)'];

% Kenar uzunluklari:
for i = 1 : size(edge, 1)
    KU.G(i,1) = abs(lat2eqdist(e, B(edge(i, 1))) - lat2eqdist(e, B(edge(i, 2))));
    if B(edge(i, 1)) == B(edge(i, 2))
        KU.Sp(i, 1) = longt2dist(e, L(edge(i, 1)), L(edge(i, 2)), B(edge(i, 1)));
        [~, ~, KU.S(i, 1)] = JTP2(e, B(edge(i, 1)), B(edge(i, 2)), L(edge(i, 1)), L(edge(i, 2)));
        KU.diff(i, 1) = KU.Sp(i, 1) - KU.S(i, 1);
    end
end

KENAR = round([double(edge) ,KU.G, KU.Sp, KU.S, KU.diff], 5);

% Pafta Alani
for i = 1 : size(area,1)
    PA.F(i, 1) = ellipsoidArea(e, B(area(i, 1)), B(area(i, 3)), L(area(i, 1)), L(area(i, 3)));
    [PA.km(i, 1), PA.hektar(i, 1), PA.donum(i, 1)] = areaConvert(PA.F(i, 1));
end

ALAN = round([double(area), PA.F, PA.km, PA.hektar, PA.donum], 5);

% Agirlik merkezi yaricaplari
for i = 1: size(area, 1)
    MA.Bort(i, 1) = (B(area(i, 1)) + B(area(i, 3))) / 2;
    [MA.N(i, 1), MA.M(i, 1), MA.R(i, 1)] = radiusCurveture(e, MA.Bort(i, 1));
end

AGIRLIK = round([double(area), MA.N, MA.M, MA.R], 3);


% Jeosantrik koordinatlar
[x, y, z] = geog2geoc(e, B(cidx), L(cidx));
JEOCENTRIC = [double(cidx)', x', y', z'];


% terminal log ciktisini yazdirma

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

% tablolari excel'e yazdirma fonksiyonu
excelWriter = main{2};

% Alanları cevirme
function [km, hektar, donum] = areaConvert(m)
    km = m / 1e6;
    hektar = m / 1e4;
    donum =  m / 1e3;
end

% tablolari excel'e yazdirma fonksiyonu
%   args:
%       content: table // AGIRLIK  etc.
%       filename: str // grs80-agirlik-50 etc.
%   returns:
%       excel file.
function writeExcel(content, filename)
    T = table(content);
    filename = [filename, '.xlsx'];
    writetable(T,filename,'Sheet','MyNewSheet','WriteVariableNames',false);
end