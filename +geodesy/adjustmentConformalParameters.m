%   args:
%       inverse: transforming coordinate system,
%       forward: transforming into,
%   returns:
%       t: (struct) conformal parameters,
function t = adjustmentConformalParameters(inverse, forward)
%     G = shiftAndNorm(inverse);

    giveMat = @(xy) [xy(1) -xy(2) 1 0;
                     xy(2)  xy(1) 0 1];
                 
    giveFwd = @(xy) xy';
                 
    L = giveFwd(forward);
    
    out.A  = [];
                 
    for n = 1 : size(inverse, 1)
        out.A = [out.A; giveMat(inverse(n, :))];
    end
    A = [out.A];
    
    N = (A' * A);
    n = A' * L(:);
    
%     E = G * G'; 
%     Np = (N + E)^-1 - E;

    x = N \ n;
    
    clc,
    
    s = sqrt(x(1)^2 + x(2)^2);
    
    beta = atan(x(2) / x(1));
    
    t.rotateMatrix = [ x(1), x(2); 
                      -x(2), x(1)];
                  
    t.shiftVector = [x(3); x(4)];
    t.scaleFactor = s;
    t.rotateAngle = beta;
end

% - HELPER FUNCTION -
% function G = shiftAndNorm(inverse)
%     square = @(x, y) sqrt(x.^2 + y.^2);
%     
%     len = size(inverse, 1);
%     
%     giveMat = @(xy) [1/sqrt(len)  0          -xy(2)   xy(1);
%                      0          1/sqrt(len)   xy(1)   xy(2)];
%                  
%     xs = mean(inverse(:, 1));
%     ys = mean(inverse(:, 2));
%     
%     x_ = inverse(:, 1) - xs;
%     y_ = inverse(:, 2) - ys;
%     
%     x__ = x_ ./ square(x_, y_);
%     y__ = y_ ./ square(x_, y_);
% 
%     % normed coordinate,
%     
%     xy = [x__, y__];
%     
%     out.G = [];
%     for n = 1 : len
%         out.G = [out.G; giveMat(xy(n, :))];
%     end
%     G = [out.G];
% end