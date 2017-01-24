classdef LoG < SynEM.Feature.TextureFeature
    %LOG Multi-dimensional laplacian of gaussian filter.
    %
    % PROPERTIES
    % sigma: [Nx1] array of double specifying the standard deviation
    %       of the gaussian kernel for each dimension. Each value must be
    %       bigger equal to zero and a zero corresponds to no filtering in
    %       the respective dimension.
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
        convMode = 'same';
    end
    
    methods
        function obj = LoG(sigma, truncate, convMode)
            obj.name = 'LoG';
            obj.sigma = sigma;
            if exist('truncate','var') && ~isempty(truncate)
                obj.truncate = truncate;
            end
            if exist('convMode','var') && ~isempty(convMode)
                obj.convMode = convMode;
            end
            obj.numChannels = 1;
            obj.border = 2.*ceil(obj.sigma*obj.truncate);
        end

        function fm = calc(obj, raw)
            fm = obj.calculate(raw, obj.sigma, obj.truncate, obj.convMode);
        end
    end
    
    methods (Static)
        function fm = calculate(raw, sigma, truncate, convMode)
            fm = zeros(size(raw),'like',raw);
            nD = ndims(raw);
            for dim = 1:nD
                order = zeros(nD,1);
                order(dim) = 2;
                fm = fm + SynEM.Feature.GaussFilter.calculate(raw, ...
                    sigma, truncate, order, convMode);
            end
        end
    end
    
end

