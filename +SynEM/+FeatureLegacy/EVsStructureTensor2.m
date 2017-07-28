classdef EVsStructureTensor2 < SynEM.Feature.TextureFeature
    %EVSSTRUCTURETENSOR Eigenvalues for structure tensor sorted by
    % increasing absolute value.
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

    properties
        sigmaW
        filterSizW
        sigmaD
        filterSizD
    end

    methods
        function obj = EVsStructureTensor2(sigmaW, filterSizW, ...
                sigmaD, filterSizD)
            obj.name = 'EVsStructureTensor';
            obj.sigmaW = sigmaW;
            obj.filterSizW = filterSizW;
            obj.sigmaD = sigmaD;
            obj.filterSizD = filterSizD;
            obj.numChannels = 3;
            obj.border = 2.*filterSizW + 2.*filterSizD;
        end

        function fm = calc(obj, raw)
            %calculate structure tensor
            sigW = num2cell(obj.sigmaW);
            coords = arrayfun(@(siz)(-siz:siz).', obj.filterSizW, ...
                'UniformOutput',false);
            coords = cellfun(@(coords,dim) ...
                permute(coords,[2:dim,1,(dim+1):3]), ...
                coords,num2cell(1:3),'UniformOutput',false);
            gauss = cellfun(@(coords,sigma)exp(-(coords./sigma).^2./2), ...
                coords,sigW,'UniformOutput',false);
            gauss = cellfun(@(gauss)gauss./sum(gauss),gauss, ...
                'UniformOutput',false);
            grad = obj.GaussGradient(raw, obj.sigmaD, obj.filterSizD);
            S = zeros(6, numel(grad{1}), 'like', raw);
            curRow = 1;
            for dim1 = 1:3
                for dim2 = dim1:3
                    tmp = grad{dim1}.*grad{dim2};
                    for dim3=1:3
                        tmp = convn(tmp, gauss{dim3}, 'same');
                    end
                    S(curRow, :) = tmp(:);
                    curRow = curRow + 1;
                end
            end
            [nx,ny,nz] = size(raw);
            ev = SynEM.Aux.eig3S(S);
            clear S

            %sort evs by absolute value
            SynEM.Aux.sortAbs(ev);
            fm = cell(3,1);
            fm{1} = reshape(ev(1,:),nx,ny,nz);
            fm{2} = reshape(ev(2,:),nx,ny,nz);
            fm{3} = reshape(ev(3,:),nx,ny,nz);
        end
    end
    methods (Static)
        function grad = GaussGradient( raw, sigma, filterSize )
        %GAUSSGRADIENTMAGNITUDE Gauss gradient of a 3D volume.
        % INPUT sigma: Array containing standard deviation used for
        %       	gaussian filterin each dimension.
        %       filterSize: Size of the resulting filter in each dimension
        %       	and direction. The resulting filter has a size of
        %       	2*filterSize + 1 and a total boundary of 2*filterSize.
        % OUTPUT grad: cell array of length 3 containing the gradient for
        %        	each dimension.
        % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

        sigma = num2cell(sigma);
        coords = arrayfun(@(siz)(-siz:siz).', filterSize, ...
            'UniformOutput',false);
        coords = cellfun(@(coords,dim)permute( ...
            coords,[2:dim,1,(dim+1):3]),coords,num2cell(1:3), ...
            'UniformOutput',false);
        gauss = cellfun(@(coords,sigma)exp(-(coords./sigma).^2./2), ...
            coords,sigma,'UniformOutput',false);
        gauss = cellfun(@(gauss)gauss./sum(gauss),gauss, ...
            'UniformOutput',false);
        gaussD = cellfun(@(coords,gauss,sigma) ...
            -gauss.*coords./(sigma.^2),coords,gauss, sigma, ...
            'UniformOutput',false);
        grad = cellfun(@(gaussD)convn(raw,gaussD,'same'), gaussD, ...
            'UniformOutput',false);
        for dim1 = 1:3
            for dim = setdiff(1:3,dim1)
                grad{dim1} = convn(grad{dim1}, gauss{dim}, 'same');
            end
        end
        end
    end

end
