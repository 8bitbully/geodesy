%   args:
%       name: ellipsoid name
%   returns:
%       ellipsoid: ReferenceEllipsoid object,
function ellipsoid = refEllipsoid(name)
    ellipsoid = src.ReferenceEllipsoid(name);
end