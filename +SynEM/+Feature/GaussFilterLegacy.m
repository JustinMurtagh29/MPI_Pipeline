classdef GaussFilterLegacy < SynEM.Feature.TextureFeature
    %GaussFilterLegacy Multidimensional gaussian filter.
    %
    % PROPERTIES
    % sigma: [3x1] array of double specifying the standard deviation
    %       of the gaussian kernel for each dimension. Each value must be bigger
    %       equal to zero and a zero corresponds to no filtering in the
    %       respective dimension.
    % filterSiz: [3x1] double
    %       Filter size of gaussian filter in each dimension.
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

    properties
        sigma
        filterSiz
    end

    methods
        function obj = GaussFilterLegacy(sigma, filterSiz)
            obj.name = 'GaussFilt';
            obj.sigma = sigma;
            obj.filterSiz = filterSiz;
            obj.border = 2.*filterSiz;
            obj.numChannels = 1;
        end

        function fm = calc(obj, raw)
            fm = obj.calculate(raw, obj.sigma, obj.filterSiz);
        end
    end

    methods (Static)
        function fm = calculate(raw, sigma, filterSize)
            sigma = num2cell(sigma);
            coords = arrayfun(@(siz)(-siz:siz).', filterSize, ...
                'UniformOutput',false);
            coords = cellfun(@(coords,dim)permute(coords, ...
                [2:dim,1,(dim+1):3]),coords,num2cell(1:3), ...
                'UniformOutput',false);
            gauss = cellfun(@(coords,sigma)exp(-(coords./sigma).^2./2), ...
                coords,sigma,'UniformOutput',false);
            gauss = cellfun(@(gauss)gauss./sum(gauss), gauss, ...
                'UniformOutput',false);
            fm = raw;
            for dim = 1:3
                fm = convn(fm,gauss{dim},'same');
            end
        end
    end
end
