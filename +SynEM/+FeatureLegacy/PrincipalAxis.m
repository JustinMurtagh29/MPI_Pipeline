classdef PrincipalAxis < SynEM.Feature.ShapeFeature
    %PRINCIPALAXIS Principal axis feature.
    %
    % PROPERTIES
    % mode: String specifying the computation mode
    %       'length': Length of the principal axis for each shape ordered
    %           in decreasing magnitude.
    %       'prod': Scalar product of first principal axis between two
    %           shapes. In this case vols needs to be of size [Nx2] and the
    %           product is calculated for the first principal axis of the
    %           shapes in each row.
    %       'eccentricity': Todo
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

    properties
        mode
        voxelSize = [];
    end

    methods
        function obj = PrincipalAxis(mode, voxelSize, numFeat)
            %INPUT mode: String specifying the computation mode
            %       'length': Get the principal axis length for all
            %           subvolumes.
            %       'product': Get the product of the main principle axis.
            %           In this case the vol input of calculate must be a
            %           [Nx2] cell array.
            %      voxelSize: (Optional) [1x3] double specifying the voxel
            %           size for each dimension.
            %           (Default: [], incoming structure is not modified)
            %      numFeat: (Optional) Number of features calculated by
            %           Volume.
            %           (Default: depending on mode - see below)

            switch mode
                case 'length'
                    if exist('numFeat','var') && ~isempty(numFeat)
                        obj.numFeatures = numFeat;
                    else
                        obj.numFeatures = 3;
                    end
                case 'prod'
                    obj.numFeatures = 1;
                case 'eccentricity'
                    error('Mode %s not yet implemented', mode);
                otherwise
                    error('Unknown mode %s', mode);
            end
            obj.mode = mode;
            obj.name = 'PrAx';
            if exist('voxelSize','var') && ~isempty(voxelSize)
                obj.voxelSize = voxelSize;
            end

        end

        function feat = calculate(obj, vol, selFeat)
            if ~isempty(obj.voxelSize)
                vol = cellfun(@(x)bsxfun(@times,x,obj.voxelSize),vol, ...
                    'UniformOutput',false);
            end
            switch obj.mode
                case 'length'
                    feat = cell2mat(cellfun(@obj.paLength,vol, ...
                        'UniformOutput',false));
                    if ~exist('selFeat','var') || isempty(selFeat)
                        selFeat = true(obj.numFeatures,1);
                    end
                    feat = feat(:,selFeat);
                case 'prod'
                    if exist('selFeat','var') && ~selFeat
                        feat = [];
                        return;
                    end
                    feat = cellfun(@obj.paProduct,vol(:,1),vol(:,2));
            end
        end

        function L = paLength(~, X)
            try
                [~,~,L] = pca(X,'Algorithm','eig');
            catch
                [~,~,L] = pca(X,'Algorithm','svd');
                L(end:3) = 0;
            end
            if iscolumn(L)
                L = L';
            end
        end
        function p = paProduct(~, X, Y)
            if size(X, 1) < 4
                X = repmat(X, ceil(4/size(X, 1)), 1);
            end
            if size(Y, 1) < 4
                X = repmat(Y, ceil(4/size(Y, 1)), 1);
            end
            try %legacy
                pca1 = pca(X,'Algorithm','eig');
            catch
                %in case of not positive semidefinite
                pca1 = pca(X,'Algorithm','svd');
            end
            try %legacy
                pca2 = pca(Y,'Algorithm','eig');
            catch
                %in case of not positive semidefinite
                pca2 = pca(Y,'Algorithm','svd');
            end
            p = pca1(:,1)'*pca2(:,1);
        end
    end

end
