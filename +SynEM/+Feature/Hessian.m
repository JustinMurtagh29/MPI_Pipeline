classdef Hessian < SynEM.Feature.TextureFeature
    %HESSIAN Multi-dimensional hessian matrix.
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
        function obj = Hessian(sigma, truncate, convMode)
            obj.name = 'Hessian';
            obj.sigma = sigma;
            if exist('truncate','var') && ~isempty(truncate)
                obj.truncate = truncate;
            end
            if exist('convMode','var') && ~isempty(convMode)
                obj.convMode = convMode;
            end
            obj.numChannels = length(sigma)*(length(sigma) - 1)/2;
            obj.border = 2.*ceil(obj.sigma*obj.truncate);
        end
        
        function H = calc(obj, raw)
            H = obj.calculate(raw, obj.sigma, obj.truncate, obj.convMode);
        end
    end
    
    methods (Static)
        function H = calculate(raw, sigma, truncate, convMode)
            nD = ndims(raw);
            H = cell(nD,nD);
            for i = 1:nD
                for j = i:nD
                    order = zeros(ndims(raw),1);
                    order(i) = order(i) + 1;
                    order(j) = order(j) + 1;
                    H{i,j} = SynEM.Feature.GaussFilter.calculate(raw, ...
                        sigma, truncate, order, convMode);
                end
            end
        end
    end
end

