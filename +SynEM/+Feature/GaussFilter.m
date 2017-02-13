classdef GaussFilter < SynEM.Feature.TextureFeature
    %GaussFilter Multidimensional gaussian filter.
    %
    % PROPERTIES
    % sigma: [1x1] or [Nx1] array of double specifying the standard deviation
    %       of the gaussian kernel for each dimension. Each value must be bigger
    %       equal to zero and a zero corresponds to no filtering in the
    %       respective dimension.
    % truncate: Double. Filter is truncated at ceil(truncate.*sigma) many
    %       standard deviations.
    %       (Default: 3)
    % order: [1x1] or [Nx1] array of integer specifying the order of the filter
    %       in the respective dimension or for al dimensions. 0 order
    %       corresponds to a regular gaussian filter and higher order correspond
    %       to the convolution with derivatives of the gaussian filter.
    %       Currently only order 1 and 2 are supported.
    %       (Default: 0)
    % convMode: String specifying the convolution mode.
    %       (see convn).
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

    properties
        sigma
        truncate = 3;
        order = 0;
        convMode = 'same';
    end

    methods
        function obj = GaussFilter(sigma, truncate, order, convMode)
            obj.name = 'GaussFilt';
            obj.sigma = sigma;
            if exist('truncate','var') && ~isempty(truncate)
                obj.truncate = truncate;
            end
            if exist('convMode','var') && ~isempty(convMode)
                obj.convMode = convMode;
            end
            if exist('order','var') && ~isempty(order)
                obj.order = order;
            end
            obj.numChannels = 1;
            obj.border = 2.*ceil(obj.sigma*obj.truncate);
        end

        function fm = calc(obj, raw)
            fm = obj.calculate(raw, obj.sigma, obj.truncate, obj.order, ...
                obj.convMode);
        end
    end

    methods (Static)
        function fm = calculate(raw, sigma, truncate, order, convMode)
            numDims = ndims(raw);
            if length(sigma) == 1
                sigma = repmat(sigma, numDims, 1);
            end
            if length(order) == 1
                order = repmat(order, numDims, 1);
            end
            fm = raw;

            %separable filter
            for dim = 1:numDims
                if sigma(dim) == 0
                    continue;
                end
                coords = ceil(sigma(dim)*truncate);
                sz = ones(1,numDims);
                sz(dim) = 2*coords + 1;
                x = reshape(-coords:coords, sz);
                h = exp(-x.^2/2./sigma(dim).^2);
                h = h./sum(h(:));
                switch order(dim)
                    case 0
                        %already calculated
                    case 1
                        h = -x./sigma(dim)^2.*h;
                    case 2
                        h = (x.^2./sigma(dim)^2 - 1)./sigma(dim).^2.*h;
                    otherwise
                        error('An order of ''%d'' is not implemented.', ...
                            order(dim));
                end
                fm = convn(fm,h,convMode);
            end
        end
    end
end
