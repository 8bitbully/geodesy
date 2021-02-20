%
% inverse: [X1, Y1;
%           X2, Y2]
%
% forward: [X1, Y1;
%           X2, Y2]

%   args:
%       inverse: transforming coordinate system,
%       forward: transforming into,
%   returns:
%       t: (struct) conformal parameters,
function t = conformalParameters(inverse, forward)
    giveMat = @(xy) [xy(1)  xy(2) 1 0;
                     xy(2) -xy(1) 0 1];
                 
    giveFwd = @(xy) xy';
                 
    rightVec = giveFwd(forward);
    
    params = inv([giveMat(inverse(1, :)); giveMat(inverse(2, :))]) * rightVec(:);
    
    s = sqrt(params(1)^2 + params(2)^2);
    
    beta = atan(params(2) / params(1));
    
    t.rotateMatrix = [ params(1), params(2); 
                      -params(2), params(1)];
                  
    t.shiftVector = [params(3); params(4)];
    t.scaleFactor = s;
    t.rotateAngle = beta;
end