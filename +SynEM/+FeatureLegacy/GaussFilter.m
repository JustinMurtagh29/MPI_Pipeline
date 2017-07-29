classdef GaussFilter < SynEM.Feature.TextureFeature
    %GaussFilter 3-dimensional gaussian filter.
    %
    % PROPERTIES
    % sigma: [3x1] array of double specifying the standard deviation
    %       of the gaussian kernel for each dimension. Each value must be bigger
    %       equal to zero and a zero corresponds to no filtering in the
    %       respective dimension.
    % filterSiz: [3x1] double
    %       Filter size of gaussian filter in each dimension.
    % n_mean: double
    %       The mean used for raw data normalization.
    % n_std: double
    %       The standard deviation used for raw data normlaization.
    %
    % NOTE n_mean and n_std are used to ensure a consistent normalization
    %      w.r.t. to uint8 raw data.
    %
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

    properties
        sigma
        filterSiz
        n_mean = 0
        n_std = 1
    end

    methods
        function obj = GaussFilter(sigma, filterSiz, n_mean, n_std)
            obj.name = 'GaussFilter';
            obj.sigma = sigma;
            obj.filterSiz = filterSiz;
            obj.border = 2.*filterSiz;
            obj.numChannels = 1;
            if exist('n_mean', 'var') && ~isempty(n_mean)
                obj.n_mean = n_mean;
            end
            if exist('n_std', 'var') && ~isempty(n_std)
                obj.n_std = n_std;
            end
        end

        function fm = calc(obj, raw)
            fm = obj.calculate(raw, obj.sigma, obj.filterSiz, ...
                obj.n_mean, obj.n_std);
        end
    end

    methods (Static)
        function fm = calculate(raw, sigma, filterSize, n_mean, n_std)
            
            % normalize
            raw = raw + n_mean/n_std;
            
            % calculate filter
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
