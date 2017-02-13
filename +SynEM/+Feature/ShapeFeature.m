classdef ShapeFeature < handle
    %SHAPEFEATURE Bass class for shape features.
    %
    % PROPERTIES
    % name: String specifying the feature name
    % numFeatures: Integer specifying the total number of output features.
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

    properties
        name
        numFeatures %supposed to be public such that it can be modified
                    %based on input
    end

    methods
        %Main feature calculation method taking the volume
        % INPUT vol: [NxM] Cell array. Each cell contains a [Nx3] integer
        %           arrays of coordinates pixel coordinates. Exact use can
        %           depend on the feature.
        %       selFeat: (Optional) Logical array of length numFeatures
        %           specifying which features need to be calculated.
        %           (Default: All features).
        % OUTPUT feat: Matrix whose exact size might depend on the feature
        %           and input.
        feat = calculate(obj, vol, selFeat)
    end
end
