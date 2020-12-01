%   Ellipsoid:
%   
%   
%
%

% References:
%       Comparison of Different Algorithms between Geocentric and Geodetic Coordinates
%       Harita Dergisi Temmuz 2011 Sayı 146
%
%       https://avesis.ktu.edu.tr//yfaruk
%
%
% for the updated version : github.com/solounextracto
% @author: 
% @date: 20202811
classdef NivoEllipsoid < matlab.mixin.CustomDisplay
    methods (Access = protected)
        function header = getHeader(this)
            className = matlab.mixin.CustomDisplay.getClassNameForHeader(this);
            newHeader = [className, 'Ellipsoid: ', this.Name] ;
            header = sprintf('%s\n\t', newHeader);
        end
        function propgrp = getPropertyGroups(~)
            Title = 'Definition Parameters';
            propTitle = {'SemimajorAxis', 'SemiminorAxis', ...
                'InverseFlattening', 'Eccentricity'};
            propgrp(1) = matlab.mixin.util.PropertyGroup(propTitle, Title);
        end
        function footer = getFooter(~)
            footer = sprintf('%s\n', 'Nivo Ellipsoid');
        end
    end
end