classdef GaussGradient < SynEM.Feature.TextureFeature
    %GAUSSGRADIENT Multi-dimensional gauss-gradient filter.
    %
    % PROPERTIES
    % sigma: [Nx1] array of double specifying the standard deviation
    %       of the gaussian kernel for each dimension. Each value must be bigger
    %       equal to zero and a zero corresponds to no filtering in the
    %       respective dimension.
    % truncate: Double. Filter is truncated at ceil(truncate.*sigma) many
    %       standard deviations.
    %       (Default: 3)
    % convMode: String specifying the convolution mode.
    %       (see convn).
    %
    % NOTE The numChannels properties is set to length(sigma). If sigma is
    %      specified as a scalar manually reset numChannels if required.
    %
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

    properties
        sigma
        truncate = 3;
        convMode = 'same';
    end

    methods
        function obj = GaussGradient(sigma, truncate, convMode)
            obj.name = 'GaussGrad';
            obj.sigma = sigma;
            if exist('truncate','var') && ~isempty(truncate)
                obj.truncate = truncate;
            end
            if exist('convMode','var') && ~isempty(convMode)
                obj.convMode = convMode;
            end
            obj.numChannels = length(sigma);
            obj.border = 2.*ceil(obj.sigma*obj.truncate);
        end

        function fm = calc(obj, raw)
            fm = obj.calculate(raw, obj.sigma, obj.truncate, obj.convMode);
        end
    end

    methods (Static)
        function fm = calculate(raw, sigma, truncate, convMode)
            numDims = ndims(raw);
            fm = cell(numDims,1);
            for dim = 1:numDims
                order = zeros(numDims,1);
                order(dim) = 1;
                fm{dim} = SynEM.Feature.GaussFilter.calculate(raw, ...
                    sigma, truncate, order, convMode);
            end
        end
    end
end
