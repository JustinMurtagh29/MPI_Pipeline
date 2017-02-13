classdef TextureFeature
    %TEXTUREFEATURE Base class for texture features/image filters for
    % multidimensional images.
    %
    % PROPERTIES
    % name: String specifying the feature name
    % border: [1xN] array of integer containing the total border of the filter
    %       in each image dimension. This always corresponds to the value
    %       such that there are no boundary effects anymore (i.e. to the border
    %       of a 'valid' convolution).
    % numChannels: Integer specifying the number of channels of the feature
    %       (i.e. the number of output feature maps).
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

    properties
        name
        border
        numChannels
    end

    methods (Abstract)
        %fm should be a cell array of length numChannels containing the
        %different feature maps
        fm = calc(obj, I);
    end
end
