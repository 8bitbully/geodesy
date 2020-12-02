clearvars; clc;
% Version: 9.9.0.1467703 (R2020b)

main = localfunctions;
[kose, kenar, alan, agirlik, jeosentrik] = getFunctions;
% elipsoid secimi
e = ReferenceEllipsoid('hayford');

% pafta secimi
pf = Pafta.PAFTA_100;

% pafta indexleri
[cidx, area, edge] = paftaEdge(pf);
[B, L] = paftaParser(pf, 39.5, 40.5);

% Pafta köse koordinatlari:
PAFTA_KOSE = kose(cidx, B, L);

% Kenar uzunluklari:
KENAR = kenar(e, edge, B, L);

% Pafta Alani
ALAN = alan(e, area, B, L);

% Agirlik merkezi yaricaplari
% return: - N - R - M -
AGIRLIK = agirlik(e, area, B);

% Jeosentrik Koordinatlar
GEOCENTRIC = jeosentrik(e, cidx, B, L);

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
excelWriter = main{1};

% tablolari excel'e yazdirma fonksiyonu
%   args:
%       content: table // AGIRLIK, KENAR  etc.
%       filename: str // grs80-agirlik-50 etc.
%   returns:
%       excel file.
%
%   Example:
%   command window:
%   excelWriter('hayford-25', PAFTA_KOSE, KENAR, ALAN, AGIRLIK, GEOCENTRIC)
function writeExcel(filename, varargin)
    n = length(varargin);
%     T = table(content);
    filename = [filename, '.xlsx'];
    for i = 1 : n
%         writetable(T,filename,'Sheet','MyNewSheet','WriteVariableNames',false);
        writetable(varargin{i}, filename,'Sheet',inputname(i+1) ,'WriteVariableNames',true);
    end
end