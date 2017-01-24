classdef LoG < SynEM.Feature.TextureFeature
    %LOGL Laplacian of Gaussian filter.
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>
    
    properties
        sigma
        filterSize
    end
    
    methods
        function obj = LoG(sigma, filterSize)
            % INPUT sigma: Array containing standard deviation used for 
            %           gaussian filter in each dimension.
            %       filterSize: Size of the resulting filter in each 
            %           dimension and direction. The resulting filter has a
            %           size of 2*filterSize + 1 and a total boundary of
            %           2*filterSize.

            obj.name = 'LoG';
            obj.sigma = sigma;
            obj.filterSize = filterSize;
            obj.numChannels = 1;
            obj.border = 2.*obj.filterSize;
        end
        
        function fm = calc(obj, raw)
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
            gaussD2 = cellfun(@(coords,gauss,sigma) ...
                gauss.*(coords.^2./sigma.^2-1)./sigma.^2, ...
                coords,gauss,sig,'UniformOutput',false);
            fm = zeros(size(raw));
            for dim1 = 1:3
                tmp = convn(raw,gaussD2{dim1},'same');
                for dim = setdiff(1:3,dim1)
                    tmp = convn(tmp,gauss{dim},'same');
                end
                fm = fm + tmp;
            end
        end
    end
    
end

