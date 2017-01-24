classdef GaussGradMagnitude < SynEM.Feature.TextureFeature
    %GAUSSGRADMAGNITUDE Gauss gradient magnitude filter.
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>
    
    properties
        sigma
        filterSize
    end
    
    methods
        
        function obj = GaussGradMagnitude(sigma, filterSize)
            obj.name = 'GaussGradMagnitude';
            obj.sigma = sigma;
            obj.filterSize = filterSize;
            obj.numChannels = 1;
            obj.border = 2.*obj.filterSize;
        end
        
        function grad = gaussGradient(obj, raw)
            %GAUSSGRADIENTMAGNITUDE Calculate the gauss gradient of a 3D volume.
            % INPUT sigma: [3x1] double
            %           Standard deviation used for gaussian filter in each
            %           dimension.
            %       filterSize: [3x1] int
            %           Filter size in each dimension and direction. The
            %           resulting filter has a size of 2*filterSize + 1
            %           and a total boundary of 2*filterSize.
            % OUTPUT grad: [3x1] cell array of size(raw) arrays
            %           Gradient image stack for each direction.
            
            sig = num2cell(obj.sigma);
            coords = arrayfun(@(siz)(-siz:siz).',obj.filterSize, ...
                'UniformOutput',false);
            coords = cellfun(@(coords,dim) ...
                permute(coords,[2:dim,1,(dim+1):3]), ...
                coords,num2cell(1:3),'UniformOutput',false);
            gauss = cellfun(@(coords,sigma)exp(-(coords./sigma).^2./2), ...
                coords,sig,'UniformOutput',false);
            gauss = cellfun(@(gauss)gauss./sum(gauss),gauss, ...
                'UniformOutput',false);
            gaussD = cellfun(@(coords,gauss,sigma) ...
                -gauss.*coords./(sigma.^2), ...
                coords,gauss,sig,'UniformOutput',false);
            grad = cellfun(@(gaussD)convn(raw,gaussD,'same'), gaussD, ...
                'UniformOutput',false);
            for dim1 = 1:3
                for dim = setdiff(1:3,dim1)
                    grad{dim1} = convn(grad{dim1},gauss{dim},'same');
                end
            end
        end
        
        function fm = calc(obj, raw)
            grad = obj.gaussGradient(raw);
            fm = sqrt(grad{1}.^2 + grad{2}.^2 + grad{3}.^2);
        end
    end
    
end

