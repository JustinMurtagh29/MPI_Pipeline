classdef DoG < SynEM.Feature.TextureFeature
    %DOG Multi-dimensional difference of gaussians filter.
    % PROPERTIES
    % sigma: [Nx1] array of double specifying the standard deviation
    %       of the gaussian kernel for each dimension. Each value must be bigger
    %       equal to zero and a zero corresponds to no filtering in the
    %       respective dimension.
    % k: (Optional) The standard deviation of the second gaussian is the
    %       standard deviation of the first gaussian multiplied by the
    %       scalar k.
    %       (Default: 1.5)
    % truncate: Double. Filter is truncated at ceil(truncate.*sigma) many
    %       standard deviations.
    %       (Default: 3)
    % convMode: String specifying the convolution mode.
    %       (see convn).
    %
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>
    
    properties
        sigma
        truncate = 3;
        k = 1.5;
        convMode = 'same';
    end
    
    methods
        function obj = DoG(sigma, k, truncate, convMode)
            obj.name = 'DoG';
            obj.sigma = sigma;
            if exist('k','var') && ~isempty(k)
                obj.k = k;
            end
            if exist('truncate','var') && ~isempty(truncate)
                obj.truncate = truncate;
            end
            if exist('convMode','var') && ~isempty(convMode)
                obj.convMode = convMode;
            end
            obj.numChannels = 1;
            obj.border = 2*ceil(obj.k*obj.sigma*obj.truncate);
        end

        function fm = calc(obj, raw)
            fm = obj.calculate(raw, obj.sigma, obj.k, obj.truncate, ...
                obj.convMode);
        end
    end
    
    methods (Static)
        function fm = calculate(raw, sigma, k, truncate, convMode)
            I1 = SynEM.Feature.GaussFilter.calculate(raw, sigma, ...
                truncate, 0, convMode);
            I2 = SynEM.Feature.GaussFilter.calculate(raw, k*sigma, ...
                truncate, 0, convMode);
            fm = I1 - I2;
        end
    end
end

