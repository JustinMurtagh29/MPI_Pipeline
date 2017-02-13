classdef DoG < SynEM.Feature.TextureFeature
    %DOGL Difference of Gaussians feature.
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>
    
    properties
        sigma
        filterSize
        k = 1.5;
    end
    
    methods
        function obj = DoG(sigma, filterSize, k)
            obj.name = 'DoG';
            obj.sigma = sigma;
            obj.filterSize = filterSize;
            if exist('k','var') && ~isempty(k)
                obj.k = k;
            end
            obj.numChannels = 1;
            obj.border = 2.*obj.filterSize;
        end

        function fm = calc(obj, raw)
            sig = num2cell(obj.sigma);
            coords = arrayfun(@(siz)(-siz:siz).',obj.filterSize, ...
                'UniformOutput',false);
            coords = cellfun(@(coords,dim)permute( ...
                coords,[2:dim,1,(dim+1):3]),coords,num2cell(1:3), ...
                'UniformOutput',false);
            gauss1 = cellfun(@(coords,sigma)exp(-(coords./sigma).^2./2),...
                coords,sig,'UniformOutput',false);
            gauss1 = cellfun(@(gauss)gauss./sum(gauss),gauss1, ...
                'UniformOutput',false);
            gauss2 = cellfun(@(coords,sigma) ...
                exp(-(coords./sigma).^2./2./obj.k^2),coords,sig, ...
                'UniformOutput',false);
            gauss2 = cellfun(@(gauss)gauss./sum(gauss),gauss2, ...
                'UniformOutput',false);
            I1 = raw;
            for dim = 1:3
                I1 = convn(I1,gauss1{dim},'same');
            end
            I2 = raw;
            for dim = 1:3
                I2 = convn(I2,gauss2{dim},'same');
            end
            fm = I1 - I2;
        end
    end
    
end

