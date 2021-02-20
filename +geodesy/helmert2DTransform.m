%   args:
%       paramters: (conformalParameters::Struct) ->t
%       inverse: transforming coordinate system,
%   returns:
%       X: Transformed coordinates, [X]
%       Y: Transformed coordinates, [Y]
function [X, Y] = helmert2DTransform(parameters, inverse)
    rotateMat = parameters.rotateMatrix;
    shift = parameters.shiftVector;
%     scaleFactor = parameters.scaleFactor;
    
    convert = struct('X', [], ...
                     'Y', []);
    
    for n = 1 : size(inverse, 1)
        XY = (rotateMat * inverse(n, :)') + shift;
        convert.X(n) = XY(1);
        convert.Y(n) = XY(2);
    end
    
    X = [convert.X]';
    Y = [convert.Y]';
end