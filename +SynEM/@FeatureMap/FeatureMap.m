classdef FeatureMap < handle
    %FEATUREMAP Feature map for interfaces classification.
    %
    % PROPERTIES
    % numFeatures: Total number of features.
    % numFeaturesSelected: Total number of currently selected features.
    % featTexture: [Nx1] cell
    %       Texture feature objects.
    % featShape: [Nx1] cell
    %       Shape feature objects.
    % featShapeSubvols: [Nx1] cell
    %       Cell array of length lenght(featShape). Each cell specifies
    %       the subvolumes which are passed to the respective shape
    %       features. Possible subvolumes are 1, 2, 3 which correspond to
    %       the surface and the two largest subvolumes.
    % featShapeDirDep: [Nx1] cell of [Mx1] int
    %       Each cell contains a [Mx1] int array which specifies how the
    %       feature vector needs to be interchanged when interchanging the
    %       subvolumes of an interface, i.e. if X is the feature vector for
    %       one direction, then X(:,featShapeDirDep{i}) should be the
    %       feature vector in the other direction. This is only required
    %       for mode 'direction'.
    %       (see also addFeature)
    % border: [1x3] int of integer
    %       Total border of the feature map for each dimension
    %       (i.e. the maximum border of all features).
    %       (see also Features.TextureFeature and calculate)
    % numSubvolumes: int
    %       Integer specifying the number of subvolumes per interface.
    % subvolsSize: [Nx1] int
    %       The size of the different subvolumes in nm. This is not used
    %       directly by the feature map but solely used for documentation.
    % quantiles: [Nx1] float 
    %       Quantiles for pooling. The N quantiles correspond to the first
    %       N pooling statistics. The last reamining pooling statistics 
    %       are specified in moments (see below).
    % moments: [4x1] bool
    %       Tthe (centralized) moments to calculate. The flags correspond
    %       to ['mean','var','skew','kurtosis'].
    % areaT: double
    %       Area threshold on the number of voxels of an interface surface.
    %       Only interfaces larger than the threshold are considered
    %       (i.e. intSize > areaT).
    % names: [Nx1] cell array string
    %       Names of each single feature in the feature representation of
    %       an interface.
    %       If empty it can be calculated via string via createNamesString
    % selectedFeat:  [Nx1] cell array
    %       Cell array of length(featTexture) + length(featShape). Each
    %       cell contains a array which specifies which features should be
    %       For texture features this corresponds to the the subvolumes and
    %       pooling statistics and should be saved as a matrix of size
    %       'numSummaryStatistics x (numSubvolumes + 1)'.
    %       For shape features this is a row vector of size
    %       featShape{i}.numFeatures
    %       This property is calculated via setSelectedFeat.
    % mode: string
    %       Calculation mode.
    %       'direction': (Default) Each interface is considered twice with
    %           interchanged subsegments.
    %       'single': Each interface is only considered once with the
    %           direction given by its input.
    % voxelSize: [1x3] double
    %       Voxel size in nm for each dimension.
    %       coordinates for shape features a multiplied by the voxel size
    %       before passing them to the respective feature.
    %       (Default: [1, 1, 1])
    % fRawNorm: (Optional) function handle
    %       This function handle is called on the raw data before feature
    %       calculation. (see FeatureMap.calculate)
    %       (Default: @(x)single(x)).
    %
    % USAGE
    %   see InterfaceClassification.makeFeatureMap for examples how to
    %       create feature map.
    %   To calculate the feature map use the calculate function.
    %
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>
    
    properties
        numFeatures = 0;
        numFeaturesSelected = 0;
        featTexture = cell(0,1);
        featShape = cell(0,1);
        featShapeSubvols = cell(0,1);
        featShapeDirDep = cell(0,1);
        border = [0, 0, 0];
        numSubvolumes
        subvolsSize
        quantiles
        moments = true(4,1);
        areaT = 0;
        names = cell(0);
        selectedFeat = [];
        mode = 'direction'
        voxelSize = [1, 1, 1];
        fRawNorm = '@(x)single(x)';
    end
    
    methods
        function obj = FeatureMap(subvolsSize, areaT, quantiles, ...
                moments, mode, voxelSize, fRawNorm)
            obj.subvolsSize = subvolsSize;
            obj.numSubvolumes = length(subvolsSize);
            if exist('areaT','var') || ~isempty(areaT)
                obj.areaT = areaT;
            end
            if exist('quantiles','var') && ~isempty(quantiles)
                obj.quantiles = quantiles;
            end
            if exist('moments','var') && ~isempty(moments)
                obj.moments = moments;
            end
            if exist('mode','var') && ~isempty(mode)
                switch mode
                    case {'single','direction'}
                        obj.mode = mode;
                    otherwise
                        error('Unknown mode %s', mode);
                end
            end
            if exist('voxelSize','var') && ~isempty(voxelSize)
                obj.voxelSize = voxelSize;
            end
            if exist('fRawNorm','var') && ~isempty(fRawNorm)
                if ischar(fRawNorm)
                    obj.fRawNorm = fRawNorm;
                else
                    obj.fRawNorm = func2str(fRawNorm);
                end
            end
        end
        
        function coords = ind2sub(obj, siz, IND)
            %Auxiliary function doing basically doing ind2sub and scaling
            %by voxel size.
            [x,y,z] = ind2sub(siz,IND);
            coords = [x,y,z];
            coords = bsxfun(@times,coords,obj.voxelSize);
        end
        
        function obj = saveobj(obj)
            %delete names for saving
            obj.names = {};
        end
    end
    
    methods (Access = private)
        function y = calcSumStatsInternal(obj, x, selFeat)
            %Calculate the selected summary statistics
            % x: [NxM] Cell array containing the sample statistics.
            % selFeatures: [Nx1] logical indicating the summary statistics
            %   that should be calculated. N = length(obj.quantiles) + 4
            %   (4 moment based summary statistics).
            
            y = zeros(1,sum(selFeat),'like',x);
            nQ = length(obj.quantiles);
            sQ = sum(selFeat(1:nQ));
            
            %calculate selected quantiles (currently linked to matlabs
            %quantiles function)
            if any(selFeat(1:nQ))
                y(1:sQ) = quantile(x, obj.quantiles(selFeat(1:nQ)));
            end
            
            %calculate selected moments
            selMoments = obj.moments;
            selMoments(selMoments) = selFeat(nQ + 1:end);
            y(sQ + 1:end) = ...
                SynEM.FeatureMap.moment(...
                x, selMoments);
        end
    end
    
    methods (Static)
        function y = quantile(x, p)
            %Fast quantile function without checks or interpolation.
            if isempty(p)
                y = [];
            else
                L = length(x);
                ind = max(floor(p.*L),1);
                y = x(ind);
            end
        end
        
        function y = moment(x, ind)
            %Function to calculate mean, variance, skewness and kurtosis
            %trying to reuse calculations.
            if ~any(ind)
                y = [];
            else
                
%                 %old implementation
%                 s(1) = mean(x);
%                 s(2) = var(x);
%                 s(3) = skewness(x);
%                 s(4) = kurtosis(x);
%                 y = s(ind);
                
                %fast implementation
                s = zeros(1,4);
                lx = length(x);
                s(1) = sum(x)/lx; %mean
                if ind(2) || ind(3) || ind(4) %variance
                    x0 = bsxfun(@minus,x,s(1));
                    s(2) = sum(x0.^2);
                    s2 = s(2) / lx;
                    s(2) = s(2)/(lx - 1);
                end
                if ind(3) && s2 > 0 %skewness (0 for delta distribution)
                    m3 = sum(x0.^3)/lx;
                    s(3) = m3/s2^(1.5);
                end
                if ind(4) && s2 > 0 %kurtosis (0 for delta distribution)
                    m4 = sum(x0.^4)/lx;
                    s(4) = m4/s2^2;
                end
                y = s(ind);
            end
        end
    end
end
