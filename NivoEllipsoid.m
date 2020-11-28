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