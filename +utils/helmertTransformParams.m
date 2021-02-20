function [params] = helmertTransformParams(inverse, forward)
    distance = @(p1, p2) sqrt((p1(1) - p2(1))^2 + (p1(2) - p2(2))^2);
    
    % scale factor
    lambda = distance(inverse(1, :), inverse(2, :)) / distance(forward(1, :), forward(2, :));
%     lambda = (distance(inverse(1, :)) - distance(inverse(2, :))) / ...
%         (distance(forward(1, :)) - distance(forward(2, :)));
    
    rotationAngle = @src.semt;
    
    % To: Angle of rotation,
    T = rotationAngle(inverse(1, :), inverse(2, :)) - ...
        rotationAngle(forward(1, :), forward(2, :));
    T = T * pi / 200;
    
    rotationMatrix = [cos(T) -sin(T); sin(T) cos(T)];
    
    % shifted measure,
    dXY = [forward(1, 1); forward(1, 2)] - ...
        (lambda * rotationMatrix * [inverse(1, 1); inverse(1, 2)]);
    
    params.rotationMatrix = rotationMatrix;
    params.scaleFactor = lambda;
    params.shift = dXY;
end