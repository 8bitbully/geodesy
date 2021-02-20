%   reference:
%       Bektas, Matematik Jeodezi, p. 218
%   args:
%       L: Longitude,
%       degrees: Slice degrees,
%   returns:
%       L0: Longitude of the central meridian,
%       DN: Zone number
function [L0, varargout] = dom(L, degrees)
    if '3' == degrees
        L0 = 3 * fix((L + 1.5) / 3);
    elseif '6' == degrees
        varargout{1} = fix(L / 6) + 31; % DN;
        L0 = 6 * fix(L / 6) + 3;
    else
       error('Slice width incompatible.') 
    end
end