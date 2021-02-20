%   args:
%       inv: P1 Coordinates, [X, Y]
%       fwd: P2 Coordinates, [Xi Y]
%   returns:
%       alpha: semt angle, (->grad)
function alpha = semt(inv, fwd)
    
    dY = (fwd(2) - inv(2));
    dX = (fwd(1) - inv(1));
    
    degree = (atan(dY / dX)) * 200 / pi;
    
    if 0 < dY && 0 < dX
        alpha = degree;
    elseif 0 < dY && 0 > dX
        alpha = 200 + degree;
    elseif 0 > dY && 0 > dY
        alpha = 200 + degree;
    else
        alpha = 400 + degree;
    end
end