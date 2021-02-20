function [X, Y] = helmertTransform2d(params, x, y)
    rotationMatrix = params.rotationMatrix;
    lambda = params.scaleFactor;
    dXY = params.shift;
    
    coordinate = struct('X', [], ...
                        'Y', []);
    
    for i = 1 : numel(x)
        fwdXY = dXY + lambda * (rotationMatrix * [x(i); y(i)]);
        coordinate.X(i) = fwdXY(1);
        coordinate.Y(i) = fwdXY(2);
    end
    X = [coordinate.X];
    Y = [coordinate.Y];
end