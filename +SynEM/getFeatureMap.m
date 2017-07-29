function fm = getFeatureMap( name, varargin )
%GETFEATUREMAP Create a predefined feature map.
% INPUT name: string
%           Name of the feature map (see below).
%       varargin: Any arguments required for specific feature maps.
%           For 'paper' varargin{1} specifies the voxel size in nm.
% OUTPUT fm: SynEM.FeatureMap object
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>


import SynEM.FeatureMap
import SynEM.Feature.*

switch name
    case 'paper'
        areaT = 150;
        subvolsSize = [40 80 160];
        quantiles = [0.25 0.5 0.75 0 1];
        if ~isempty(varargin)
            voxelSize = varargin{1};
        else
            voxelSize = [11.24, 11.24, 28];
        end
        moments = true(4,1);
        sigma = 12./voxelSize;
        fS = ceil(2.*sigma);
        fRawNorm = '@(x)single(x)';

        %construct feature map
        fm = FeatureMap(subvolsSize, areaT, ...
            quantiles, moments, 'direction', [1, 1, 1], fRawNorm);

        %add features
        fm.addFeature(SynEM.FeatureLegacy.Identity());
        fm.addFeature(SynEM.FeatureLegacy.EVsStructureTensor(sigma, fS, ...
            sigma, fS));
        fm.addFeature(SynEM.FeatureLegacy.EVsStructureTensor(sigma, fS, ...
            2.*sigma, 2.*fS));
        fm.addFeature(SynEM.FeatureLegacy.EVsStructureTensor(2.*sigma, ...
            2.*fS, sigma, fS));
        fm.addFeature(SynEM.FeatureLegacy.EVsStructureTensor(2.*sigma, ...
            2.*fS, 2.*sigma, 2.*fS));
        fm.addFeature(SynEM.FeatureLegacy.EVsStructureTensor(3.*sigma, ...
            3.*fS, 3.*sigma, 3.*fS));
        fm.addFeature(SynEM.FeatureLegacy.EVsHessian(sigma, fS));
        fm.addFeature(SynEM.FeatureLegacy.EVsHessian(2.*sigma, 2.*fS));
        fm.addFeature(SynEM.FeatureLegacy.EVsHessian(3.*sigma, 3.*fS));
        fm.addFeature(SynEM.FeatureLegacy.EVsHessian(4.*sigma, 4.*fS));
        fm.addFeature(SynEM.FeatureLegacy.GaussFilter(sigma, fS));
        fm.addFeature(SynEM.FeatureLegacy.GaussFilter(2.*sigma, 2.*fS));
        fm.addFeature(SynEM.FeatureLegacy.GaussFilter(3.*sigma, 3.*fS));
        fm.addFeature(SynEM.FeatureLegacy.DoG(sigma, fS, 1.5));
        fm.addFeature(SynEM.FeatureLegacy.DoG(sigma, fS, 2));
        fm.addFeature(SynEM.FeatureLegacy.DoG(2.*sigma, 2.*fS, 1.5));
        fm.addFeature(SynEM.FeatureLegacy.DoG(2.*sigma, 2.*fS, 2));
        fm.addFeature(SynEM.FeatureLegacy.DoG(3.*sigma, 3.*fS, 1.5));
        fm.addFeature(SynEM.FeatureLegacy.LoG(sigma, fS));
        fm.addFeature(SynEM.FeatureLegacy.LoG(2.*sigma, 2.*fS));
        fm.addFeature(SynEM.FeatureLegacy.LoG(3.*sigma, 3.*fS));
        fm.addFeature(SynEM.FeatureLegacy.LoG(4.*sigma, 4.*fS));
        fm.addFeature(SynEM.FeatureLegacy.GaussGradMagnitude(sigma, fS));
        fm.addFeature(SynEM.FeatureLegacy.GaussGradMagnitude(2.*sigma, 2.*fS));
        fm.addFeature(SynEM.FeatureLegacy.GaussGradMagnitude(3.*sigma, 3.*fS));
        fm.addFeature(SynEM.FeatureLegacy.GaussGradMagnitude(4.*sigma, 4.*fS));
        fm.addFeature(SynEM.FeatureLegacy.GaussGradMagnitude(5.*sigma, 5.*fS));
        fm.addFeature(SynEM.FeatureLegacy.StdFilter([5, 5, 5]));
        fm.addFeature(SynEM.FeatureLegacy.EntropyFilter([5, 5, 5]));
        fm.addFeature(SynEM.FeatureLegacy.IntVar([3, 3, 3]));
        fm.addFeature(SynEM.FeatureLegacy.IntVar([5, 5, 5]));
        fm.addFeature(SynEM.FeatureLegacy.AverageFilter('ball',3,[1, 1, 2]));
        fm.addFeature(SynEM.FeatureLegacy.AverageFilter('ball',6,[1, 1, 2]));
        fm.addFeature(SynEM.FeatureLegacy.Volume([], 3),[1, 2, 3],[1 3 2]);
        fm.addFeature(SynEM.FeatureLegacy.Volume( ...
            '@(x)single(nthroot(6*x*(11.24*11.24*28)*pi,3))',1), ...
            1, 1); %diameter from old fm
        fm.addFeature(SynEM.FeatureLegacy.PrincipalAxis('length', ...
            [1, 1, 28/11.24]), 1, [1, 2, 3]);
        fm.addFeature(SynEM.FeatureLegacy.PrincipalAxis('prod', ...
            [1, 1, 28/11.24]), [2 3], 1);
        fm.addFeature(SynEM.FeatureLegacy.ConvexHull(3),[1, 2, 3],[1 3 2]);
        
    case 'paper_opt'
        % use the "optimized" features
        areaT = 150;
        subvolsSize = [40 80 160];
        quantiles = [0.25 0.5 0.75 0 1];
        if ~isempty(varargin)
            voxelSize = varargin{1};
        else
            voxelSize = [11.24, 11.24, 28];
        end
        moments = true(4,1);
        sigma = 12./voxelSize;
        fS = ceil(2.*sigma);
        fRawNorm = '@(x)single(x)';

        %construct feature map
        fm = FeatureMap(subvolsSize, areaT, ...
            quantiles, moments, 'direction', [1, 1, 1], fRawNorm);

        %add features
        fm.addFeature(SynEM.FeatureLegacy.Identity());
        fm.addFeature(SynEM.FeatureLegacy.EVsStructureTensor2(sigma, fS, ...
            sigma, fS));
        fm.addFeature(SynEM.FeatureLegacy.EVsStructureTensor2(sigma, fS, ...
            2.*sigma, 2.*fS));
        fm.addFeature(SynEM.FeatureLegacy.EVsStructureTensor2(2.*sigma, ...
            2.*fS, sigma, fS));
        fm.addFeature(SynEM.FeatureLegacy.EVsStructureTensor2(2.*sigma, ...
            2.*fS, 2.*sigma, 2.*fS));
        fm.addFeature(SynEM.FeatureLegacy.EVsStructureTensor2(3.*sigma, ...
            3.*fS, 3.*sigma, 3.*fS));
        fm.addFeature(SynEM.FeatureLegacy.EVsHessian2(sigma, fS));
        fm.addFeature(SynEM.FeatureLegacy.EVsHessian2(2.*sigma, 2.*fS));
        fm.addFeature(SynEM.FeatureLegacy.EVsHessian2(3.*sigma, 3.*fS));
        fm.addFeature(SynEM.FeatureLegacy.EVsHessian2(4.*sigma, 4.*fS));
        fm.addFeature(SynEM.FeatureLegacy.GaussFilter(sigma, fS));
        fm.addFeature(SynEM.FeatureLegacy.GaussFilter(2.*sigma, 2.*fS));
        fm.addFeature(SynEM.FeatureLegacy.GaussFilter(3.*sigma, 3.*fS));
        fm.addFeature(SynEM.FeatureLegacy.DoG(sigma, fS, 1.5));
        fm.addFeature(SynEM.FeatureLegacy.DoG(sigma, fS, 2));
        fm.addFeature(SynEM.FeatureLegacy.DoG(2.*sigma, 2.*fS, 1.5));
        fm.addFeature(SynEM.FeatureLegacy.DoG(2.*sigma, 2.*fS, 2));
        fm.addFeature(SynEM.FeatureLegacy.DoG(3.*sigma, 3.*fS, 1.5));
        fm.addFeature(SynEM.FeatureLegacy.LoG(sigma, fS));
        fm.addFeature(SynEM.FeatureLegacy.LoG(2.*sigma, 2.*fS));
        fm.addFeature(SynEM.FeatureLegacy.LoG(3.*sigma, 3.*fS));
        fm.addFeature(SynEM.FeatureLegacy.LoG(4.*sigma, 4.*fS));
        fm.addFeature(SynEM.FeatureLegacy.GaussGradMagnitude(sigma, fS));
        fm.addFeature(SynEM.FeatureLegacy.GaussGradMagnitude(2.*sigma, 2.*fS));
        fm.addFeature(SynEM.FeatureLegacy.GaussGradMagnitude(3.*sigma, 3.*fS));
        fm.addFeature(SynEM.FeatureLegacy.GaussGradMagnitude(4.*sigma, 4.*fS));
        fm.addFeature(SynEM.FeatureLegacy.GaussGradMagnitude(5.*sigma, 5.*fS));
        fm.addFeature(SynEM.FeatureLegacy.StdFilter2([5, 5, 5]));
        fm.addFeature(SynEM.FeatureLegacy.EntropyFilter([5, 5, 5]));
        fm.addFeature(SynEM.FeatureLegacy.IntVar2([3, 3, 3]));
        fm.addFeature(SynEM.FeatureLegacy.IntVar2([5, 5, 5]));
        fm.addFeature(SynEM.FeatureLegacy.AverageFilter('ball',3,[1, 1, 2]));
        fm.addFeature(SynEM.FeatureLegacy.AverageFilter('ball',6,[1, 1, 2]));
        fm.addFeature(SynEM.FeatureLegacy.Volume([], 3),[1, 2, 3],[1 3 2]);
        fm.addFeature(SynEM.FeatureLegacy.Volume( ...
            '@(x)single(nthroot(6*x*(11.24*11.24*28)*pi,3))',1), ...
            1, 1); %diameter from old fm
        fm.addFeature(SynEM.FeatureLegacy.PrincipalAxis('length', ...
            [1, 1, 28/11.24]), 1, [1, 2, 3]);
        fm.addFeature(SynEM.FeatureLegacy.PrincipalAxis('prod', ...
            [1, 1, 28/11.24]), [2 3], 1);
        fm.addFeature(SynEM.FeatureLegacy.ConvexHull(3),[1, 2, 3],[1 3 2]);
        
    case 'paperNorm'
        areaT = 150;
        subvolsSize = [40 80 160];
        quantiles = [0.25 0.5 0.75 0 1];
        voxelSize = [11.24, 11.24, 28];
        moments = true(4,1);
        sigma = 12./voxelSize;
        fS = ceil(2.*sigma);
        n_mean = 122;
        n_std = 22;
        fRawNorm = func2str(@(x)(single(x) - 122)./22);

        %construct feature map
        fm = FeatureMap(subvolsSize, areaT, ...
            quantiles, moments, 'direction', [1, 1, 1], fRawNorm);

        %add features
        fm.addFeature(SynEM.FeatureLegacy.Identity());
        fm.addFeature(SynEM.FeatureLegacy.EVsStructureTensor(sigma, fS, ...
            sigma, fS));
        fm.addFeature(SynEM.FeatureLegacy.EVsStructureTensor(sigma, fS, ...
            2.*sigma, 2.*fS));
        fm.addFeature(SynEM.FeatureLegacy.EVsStructureTensor(2.*sigma, ...
            2.*fS, sigma, fS));
        fm.addFeature(SynEM.FeatureLegacy.EVsStructureTensor(2.*sigma, ...
            2.*fS, 2.*sigma, 2.*fS));
        fm.addFeature(SynEM.FeatureLegacy.EVsStructureTensor(3.*sigma, ...
            3.*fS, 3.*sigma, 3.*fS));
        fm.addFeature(SynEM.FeatureLegacy.EVsHessian(sigma, fS, n_mean, n_std));
        fm.addFeature(SynEM.FeatureLegacy.EVsHessian(2.*sigma, 2.*fS, n_mean, n_std));
        fm.addFeature(SynEM.FeatureLegacy.EVsHessian(3.*sigma, 3.*fS, n_mean, n_std));
        fm.addFeature(SynEM.FeatureLegacy.EVsHessian(4.*sigma, 4.*fS, n_mean, n_std));
        fm.addFeature(SynEM.FeatureLegacy.GaussFilter(sigma, fS, n_mean, n_std));
        fm.addFeature(SynEM.FeatureLegacy.GaussFilter(2.*sigma, 2.*fS, n_mean, n_std));
        fm.addFeature(SynEM.FeatureLegacy.GaussFilter(3.*sigma, 3.*fS, n_mean, n_std));
        fm.addFeature(SynEM.FeatureLegacy.DoG(sigma, fS, 1.5));
        fm.addFeature(SynEM.FeatureLegacy.DoG(sigma, fS, 2));
        fm.addFeature(SynEM.FeatureLegacy.DoG(2.*sigma, 2.*fS, 1.5));
        fm.addFeature(SynEM.FeatureLegacy.DoG(2.*sigma, 2.*fS, 2));
        fm.addFeature(SynEM.FeatureLegacy.DoG(3.*sigma, 3.*fS, 1.5));
        fm.addFeature(SynEM.FeatureLegacy.LoG(sigma, fS));
        fm.addFeature(SynEM.FeatureLegacy.LoG(2.*sigma, 2.*fS));
        fm.addFeature(SynEM.FeatureLegacy.LoG(3.*sigma, 3.*fS));
        fm.addFeature(SynEM.FeatureLegacy.LoG(4.*sigma, 4.*fS));
        fm.addFeature(SynEM.FeatureLegacy.GaussGradMagnitude(sigma, fS));
        fm.addFeature(SynEM.FeatureLegacy.GaussGradMagnitude(2.*sigma, 2.*fS));
        fm.addFeature(SynEM.FeatureLegacy.GaussGradMagnitude(3.*sigma, 3.*fS));
        fm.addFeature(SynEM.FeatureLegacy.GaussGradMagnitude(4.*sigma, 4.*fS));
        fm.addFeature(SynEM.FeatureLegacy.GaussGradMagnitude(5.*sigma, 5.*fS));
        fm.addFeature(SynEM.FeatureLegacy.StdFilter([5, 5, 5]));
        fm.addFeature(SynEM.FeatureLegacy.EntropyFilter([5, 5, 5], [], true));
        fm.addFeature(SynEM.FeatureLegacy.IntVar([3, 3, 3], n_mean, n_std));
        fm.addFeature(SynEM.FeatureLegacy.IntVar([5, 5, 5], n_mean, n_std));
        fm.addFeature(SynEM.FeatureLegacy.AverageFilter('ball',3,[1, 1, 2]));
        fm.featTexture{end}.n_mean = n_mean;
        fm.featTexture{end}.n_std = n_std;
        fm.addFeature(SynEM.FeatureLegacy.AverageFilter('ball',6,[1, 1, 2]));
        fm.featTexture{end}.n_mean = n_mean;
        fm.featTexture{end}.n_std = n_std;
        fm.addFeature(SynEM.FeatureLegacy.Volume([], 3),[1, 2, 3],[1 3 2]);
        fm.addFeature(SynEM.FeatureLegacy.Volume( ...
            '@(x)single(nthroot(6*x*(11.24*11.24*28)*pi,3))',1), ...
            1, 1); %diameter from old fm
        fm.addFeature(SynEM.FeatureLegacy.PrincipalAxis('length', ...
            [1, 1, 28/11.24]), 1, [1, 2, 3]);
        fm.addFeature(SynEM.FeatureLegacy.PrincipalAxis('prod', ...
            [1, 1, 28/11.24]), [2 3], 1);
        fm.addFeature(SynEM.FeatureLegacy.ConvexHull(3),[1, 2, 3],[1 3 2]);

    case 'paperNew'
        areaT = 150;
        subvolsSize = [40 80 160];
        quantiles = [0.25 0.5 0.75 0 1];
        voxelSize = [11.24, 11.24, 28];
        moments = true(4,1);
        sigma = 12./voxelSize;
        fRawNorm = @(x)(single(x) - 122)./22;

        %construct feature map
        fm = FeatureMap(subvolsSize, areaT, ...
            quantiles, moments, 'direction', voxelSize, fRawNorm);

        %add features
        fm.addFeature(Identity());
        fm.addFeature(EigenvaluesStructureTensor(sigma, sigma, 2, 2));
        fm.addFeature(EigenvaluesStructureTensor(sigma, 2.*sigma, 2, 2));
        fm.addFeature(EigenvaluesStructureTensor(2.*sigma, sigma, 2, 2));
        fm.addFeature(EigenvaluesStructureTensor(2.*sigma, 2.*sigma, 2, 2));
        fm.addFeature(EigenvaluesStructureTensor(3.*sigma, 3.*sigma, 2, 2));
        fm.addFeature(EigenvaluesHessian(sigma, 2));
        fm.addFeature(EigenvaluesHessian(2.*sigma, 2));
        fm.addFeature(EigenvaluesHessian(3.*sigma, 2));
        fm.addFeature(EigenvaluesHessian(4.*sigma, 2));
        fm.addFeature(GaussFilter(sigma, 2));
        fm.addFeature(GaussFilter(2.*sigma, 2));
        fm.addFeature(GaussFilter(3.*sigma, 2));
        fm.addFeature(DoG(sigma, 1.5, 2));
        fm.addFeature(DoG(sigma, 2, 2));
        fm.addFeature(DoG(2.*sigma, 1.5, 2));
        fm.addFeature(DoG(2.*sigma, 2, 2));
        fm.addFeature(DoG(3.*sigma, 1.5, 2));
        fm.addFeature(LoG(sigma, 2));
        fm.addFeature(LoG(2.*sigma, 2));
        fm.addFeature(LoG(3.*sigma, 2));
        fm.addFeature(LoG(4.*sigma, 2));
        fm.addFeature(GaussGradientMagnitude(sigma, 2));
        fm.addFeature(GaussGradientMagnitude(2.*sigma, 2));
        fm.addFeature(GaussGradientMagnitude(3.*sigma, 2));
        fm.addFeature(GaussGradientMagnitude(4.*sigma, 2));
        fm.addFeature(GaussGradientMagnitude(5.*sigma, 2));
        fm.addFeature(StdFilter([5, 5, 3]));
        fm.addFeature(EntropyFilter([5, 5, 3]));
        fm.addFeature(IntVar([3, 3, 1]));
        fm.addFeature(IntVar([5, 5, 3]));
        fm.addFeature(AverageFilter('ball', 3, voxelSize./min(voxelSize)));
        fm.addFeature(AverageFilter('ball', 6, voxelSize./min(voxelSize)));
        fm.addFeature(Volume([], 3), [1, 2, 3], [1 3 2]);
        fm.addFeature(Volume( ...
            '@(x)single(nthroot(6*x*(11.24*11.24*28)*pi,3))',1), ...
            1, 1); %diameter from old fm
        fm.addFeature(PrincipalAxis('length'), 1 ,[1, 2, 3]);
        fm.addFeature(PrincipalAxis('prod'), [2 3], 1);
        fm.addFeature(ConvexHull(3), [1, 2, 3], [1 3 2]);

    case 'revised'
        areaT = 150;
        subvolsSize = [40 80 160];
        quantiles = [0.25 0.5 0.75 0 1];
        moments = true(4,1);
        voxelSize = [11.24, 11.24, 28];
        sigma = 12./voxelSize;
        fRawNorm = @(x)(single(x) - 122)./22;
        fm = FeatureMap(subvolsSize, areaT, ...
            quantiles, moments, 'direction', voxelSize, fRawNorm);

        %add features
        fm.addFeature(Identity());
        fm.addFeature(EigenvaluesStructureTensor(sigma, sigma, 2, 2));
        fm.addFeature(EigenvaluesStructureTensor(sigma, 2.*sigma, 2, 2));
        fm.addFeature(EigenvaluesStructureTensor(2.*sigma, sigma, 2, 2));
        fm.addFeature(EigenvaluesStructureTensor(2.*sigma, 2.*sigma, 2, 2));
        fm.addFeature(EigenvaluesStructureTensor(sigma, 3.*sigma, 2, 2));
        fm.addFeature(EigenvaluesStructureTensor(2.*sigma, 3.*sigma, 2, 2));
        fm.addFeature(EigenvaluesHessian(sigma,2));
        fm.addFeature(EigenvaluesHessian(2.*sigma,2));
        fm.addFeature(EigenvaluesHessian(3.*sigma,2));
        fm.addFeature(EigenvaluesHessian(4.*sigma,2));
        fm.addFeature(GaussFilter(sigma, 2));
        fm.addFeature(GaussFilter(2.*sigma, 2));
        fm.addFeature(GaussFilter(3.*sigma, 2));
        fm.addFeature(DoG(sigma, 1.5, 2));
        fm.addFeature(DoG(sigma, 3, 2));
        fm.addFeature(DoG(sigma, 6, 2));
        fm.addFeature(DoG(sigma, 9, 2));
        fm.addFeature(DoG(2.*sigma, 1/3, 2));
        fm.addFeature(DoG(2.*sigma, 1.5, 2));
        fm.addFeature(DoG(2.*sigma, 3, 2));
        fm.addFeature(DoG(3.*sigma, 1/3, 2));
        fm.addFeature(DoG(3.*sigma, 1.5, 2));
        fm.addFeature(DoG(3.*sigma, 3, 2));
        fm.addFeature(DoG(4.*sigma, 1/4, 2));
        fm.addFeature(DoG(5.*sigma, 1/5, 2));
        fm.addFeature(LoG(sigma, 2));
        fm.addFeature(LoG(2.*sigma, 2));
        fm.addFeature(LoG(3.*sigma, 2));
        fm.addFeature(LoG(4.*sigma, 2));
        fm.addFeature(GaussGradientMagnitude(sigma, 2));
        fm.addFeature(GaussGradientMagnitude(2.*sigma, 2));
        fm.addFeature(GaussGradientMagnitude(3.*sigma, 2));
        fm.addFeature(StdFilter([5, 5, 3]));
        fm.addFeature(EntropyFilter([5, 5, 3])); %the missing normalization is intended
        fm.addFeature(EntropyFilter([5, 5, 5], [], true));
        fm.addFeature(IntVar([3, 3, 3]));
        fm.addFeature(IntVar([5, 5, 3]));
        fm.addFeature(AverageFilter('ball',3,[1, 1, 2.5]));
        fm.addFeature(AverageFilter('ball',6,[1, 1, 2.5]));
        fm.addFeature(MaximumFilter(3, sigma));
        fm.addFeature(MaximumFilter(5, sigma));
        fm.addFeature(MinimumFilter(3, sigma));
        fm.addFeature(Volume('@log', 3),[1, 2, 3],[1 3 2]);
        fm.addFeature(PrincipalAxis('length',[1, 1, 1]),1,[1, 2, 3]);
        fm.addFeature(PrincipalAxis('prod',[1, 1, 1]),[2 3],1);
        fm.addFeature(ConvexHull(3),[1, 2, 3],[1 3 2]);

    case 'small'
        areaT = 150;
        subvolsSize = [80 160];
        quantiles = [0.05 0.25 0.5 0.75 0.95];
        moments = true(4,1);
        voxelSize = [11.24, 11.24, 28];
        sigma = 12./voxelSize;
        fRawNorm = @(x)(single(x) - 122)./22;

        %construct feature map
        fm = FeatureMap(subvolsSize, areaT, ...
            quantiles, moments, 'direction', voxelSize, fRawNorm);

        %add features
        fm.addFeature(Identity());
        fm.addFeature(EigenvaluesStructureTensor(sigma, sigma, 2, 2));
        fm.addFeature(EigenvaluesStructureTensor(sigma, 2.*sigma, 2, 2));
        fm.addFeature(EigenvaluesStructureTensor(2.*sigma, sigma, 2, 2));
        fm.addFeature(EigenvaluesHessian(sigma, 2));
        fm.addFeature(EigenvaluesHessian(2.*sigma, 2));
        fm.addFeature(EigenvaluesHessian(3.*sigma, 2));
        fm.addFeature(GaussFilter(sigma, 2));
        fm.addFeature(GaussFilter(2.*sigma, 2));
        fm.addFeature(GaussFilter(3.*sigma, 2));
        fm.addFeature(DoG(sigma, 1.5, 2));
        fm.addFeature(DoG(sigma, 2, 2));
        fm.addFeature(DoG(2.*sigma, 1.5, 2));
        fm.addFeature(DoG(2.*sigma, 2, 2));
        fm.addFeature(DoG(3.*sigma, 1.5, 2));
        fm.addFeature(DoG(3.*sigma, 1/3, 2));
        fm.addFeature(LoG(sigma, 2));
        fm.addFeature(LoG(2.*sigma, 2));
        fm.addFeature(LoG(3.*sigma, 2));
        fm.addFeature(GaussGradientMagnitude(sigma, 2));
        fm.addFeature(GaussGradientMagnitude(2.*sigma, 2));
        fm.addFeature(GaussGradientMagnitude(3.*sigma, 2));
        fm.addFeature(GaussGradientMagnitude(4.*sigma, 2));
        fm.addFeature(StdFilter([5, 5, 5]));
        fm.addFeature(EntropyFilter([5, 5, 5]));
        fm.addFeature(IntVar([3, 3, 3]));
        fm.addFeature(IntVar([5, 5, 5]));
        fm.addFeature(AverageFilter('ball', 3, [1, 1, 2]));
        fm.addFeature(AverageFilter('ball', 6, [1, 1, 2]));
        fm.addFeature(Volume([], 3), [1, 2, 3], [1 3 2]);
        fm.addFeature(Volume( ...
            '@(x)single(nthroot(6*x*(11.24*11.24*28)*pi,3))',1), ...
            1, 1); %diameter from old fm
        fm.addFeature(PrincipalAxis('length'), 1 ,[1, 2, 3]);
        fm.addFeature(PrincipalAxis('prod'), [2 3], 1);
        fm.addFeature(ConvexHull(3), [1, 2, 3], [1 3 2]);

    case 'singleSubvol'
        areaT = 150;
        subvolsSize = 160;
        quantiles = [0.05 0.25 0.5 0.75 0.95];
        moments = true(4,1);
        voxelSize = [11.24, 11.24, 28];
        sigma = 12./voxelSize;
        fRawNorm = @(x)(single(x) - 122)./22;

        %construct feature map
        fm = FeatureMap(subvolsSize, areaT, ...
            quantiles, moments, 'direction', voxelSize, fRawNorm);

        %add features
        fm.addFeature(Identity());
        fm.addFeature(EigenvaluesStructureTensor(sigma, sigma, 2, 2));
        fm.addFeature(EigenvaluesStructureTensor(sigma, 2.*sigma, 2, 2));
        fm.addFeature(EigenvaluesStructureTensor(2.*sigma, sigma, 2, 2));
        fm.addFeature(EigenvaluesHessian(sigma, 2));
        fm.addFeature(EigenvaluesHessian(2.*sigma, 2));
        fm.addFeature(EigenvaluesHessian(3.*sigma, 2));
        fm.addFeature(GaussFilter(sigma, 2));
        fm.addFeature(GaussFilter(2.*sigma, 2));
        fm.addFeature(GaussFilter(3.*sigma, 2));
        fm.addFeature(DoG(sigma, 1.5, 2));
        fm.addFeature(DoG(sigma, 2, 2));
        fm.addFeature(DoG(2.*sigma, 1.5, 2));
        fm.addFeature(DoG(2.*sigma, 2, 2));
        fm.addFeature(DoG(3.*sigma, 1.5, 2));
        fm.addFeature(DoG(3.*sigma, 1/3, 2));
        fm.addFeature(LoG(sigma, 2));
        fm.addFeature(LoG(2.*sigma, 2));
        fm.addFeature(LoG(3.*sigma, 2));
        fm.addFeature(GaussGradientMagnitude(sigma, 2));
        fm.addFeature(GaussGradientMagnitude(2.*sigma, 2));
        fm.addFeature(GaussGradientMagnitude(3.*sigma, 2));
        fm.addFeature(GaussGradientMagnitude(4.*sigma, 2));
        fm.addFeature(StdFilter([5, 5, 5]));
        fm.addFeature(EntropyFilter([5, 5, 5]));
        fm.addFeature(IntVar([3, 3, 3]));
        fm.addFeature(IntVar([5, 5, 5]));
        fm.addFeature(AverageFilter('ball', 3, [1, 1, 2]));
        fm.addFeature(AverageFilter('ball', 6, [1, 1, 2]));
        fm.addFeature(Volume([], 3), [1, 2, 3], [1 3 2]);
        fm.addFeature(Volume( ...
            '@(x)single(nthroot(6*x*(11.24*11.24*28)*pi,3))',1), ...
            1, 1); %diameter from old fm
        fm.addFeature(PrincipalAxis('length'), 1 ,[1, 2, 3]);
        fm.addFeature(PrincipalAxis('prod'), [2 3], 1);
        fm.addFeature(ConvexHull(3), [1, 2, 3], [1 3 2]);
        
    case 'PaperAndCNN'
        
        %check for cnet
        if isempty(varargin) || ~isa(varargin{1}, 'Codat.CNN.cnn')
            error(['Missing input: Specify a Codat.CNN object in ' ...
                'varargin{1}.'])
        end
        cnet = varargin{1};
        
        %varargin{2} specifies which fms of the cnn to use (default last
        %hidden layer and output layer)
        if length(varargin) == 2
            fmL = varargin{2};
        else
            fmL = cnet.layer;
        end
        
        areaT = 150;
        subvolsSize = [40 80 160];
        quantiles = [0.25 0.5 0.75 0 1];
        voxelSize = [11.24, 11.24, 28];
        moments = true(4,1);
        sigma = 12./voxelSize;
        fS = ceil(2.*sigma);
        fRawNorm = '@(x)single(x)';

        %construct feature map
        fm = FeatureMap(subvolsSize, areaT, ...
            quantiles, moments, 'direction', [1, 1, 1], fRawNorm);

        %add features
        fm.addFeature(SynEM.FeatureLegacy.Identity());
        fm.addFeature(CNNFeat(cnet, ...
            struct('target_size', [64, 64, 64], 'convAlg', 'fft2', ...
            'fmL', fmL, 'display', 10), ...
            @(x)(single(x) - 122)./22));
        fm.addFeature(SynEM.FeatureLegacy.EVsStructureTensor2(sigma, fS, ...
            sigma, fS));
        fm.addFeature(SynEM.FeatureLegacy.EVsStructureTensor2(sigma, ...
            fS, 2.*sigma, 2.*fS));
        fm.addFeature(SynEM.FeatureLegacy.EVsStructureTensor2(2.*sigma, ...
            2.*fS, sigma, fS));
        fm.addFeature(SynEM.FeatureLegacy.EVsStructureTensor2(2.*sigma, ...
            2.*fS, 2.*sigma, 2.*fS));
        fm.addFeature(SynEM.FeatureLegacy.EVsStructureTensor2(3.*sigma, ...
            3.*fS, 3.*sigma, 3.*fS));
        fm.addFeature(SynEM.FeatureLegacy.EVsHessian2(sigma, fS));
        fm.addFeature(SynEM.FeatureLegacy.EVsHessian2(2.*sigma, 2.*fS));
        fm.addFeature(SynEM.FeatureLegacy.EVsHessian2(3.*sigma, 3.*fS));
        fm.addFeature(SynEM.FeatureLegacy.EVsHessian2(4.*sigma, 4.*fS));
        fm.addFeature(SynEM.FeatureLegacy.GaussFilter(sigma, fS));
        fm.addFeature(SynEM.FeatureLegacy.GaussFilter(2.*sigma, 2.*fS));
        fm.addFeature(SynEM.FeatureLegacy.GaussFilter(3.*sigma, 3.*fS));
        fm.addFeature(SynEM.FeatureLegacy.DoG(sigma, fS, 1.5));
        fm.addFeature(SynEM.FeatureLegacy.DoG(sigma, fS, 2));
        fm.addFeature(SynEM.FeatureLegacy.DoG(2.*sigma, 2.*fS, 1.5));
        fm.addFeature(SynEM.FeatureLegacy.DoG(2.*sigma, 2.*fS, 2));
        fm.addFeature(SynEM.FeatureLegacy.DoG(3.*sigma, 3.*fS, 1.5));
        fm.addFeature(SynEM.FeatureLegacy.LoG(sigma, fS));
        fm.addFeature(SynEM.FeatureLegacy.LoG(2.*sigma, 2.*fS));
        fm.addFeature(SynEM.FeatureLegacy.LoG(3.*sigma, 3.*fS));
        fm.addFeature(SynEM.FeatureLegacy.LoG(4.*sigma, 4.*fS));
        fm.addFeature(SynEM.FeatureLegacy.GaussGradMagnitude(sigma, fS));
        fm.addFeature(SynEM.FeatureLegacy.GaussGradMagnitude(2.*sigma, 2.*fS));
        fm.addFeature(SynEM.FeatureLegacy.GaussGradMagnitude(3.*sigma, 3.*fS));
        fm.addFeature(SynEM.FeatureLegacy.GaussGradMagnitude(4.*sigma, 4.*fS));
        fm.addFeature(SynEM.FeatureLegacy.GaussGradMagnitude(5.*sigma, 5.*fS));
        fm.addFeature(SynEM.FeatureLegacy.StdFilter([5, 5, 5]));
        fm.addFeature(SynEM.FeatureLegacy.EntropyFilter([5, 5, 5]));
        fm.addFeature(SynEM.FeatureLegacy.IntVar([3, 3, 3]));
        fm.addFeature(SynEM.FeatureLegacy.IntVar([5, 5, 5]));
        fm.addFeature(SynEM.FeatureLegacy.AverageFilter('ball',3,[1, 1, 2]));
        fm.addFeature(SynEM.FeatureLegacy.AverageFilter('ball',6,[1, 1, 2]));
        fm.addFeature(SynEM.FeatureLegacy.Volume([], 3),[1, 2, 3],[1 3 2]);
        fm.addFeature(SynEM.FeatureLegacy.Volume( ...
            '@(x)single(nthroot(6*x*(11.24*11.24*28)*pi,3))',1), ...
            1, 1); %diameter from old fm
        fm.addFeature(SynEM.FeatureLegacy.PrincipalAxis('length', ...
            [1, 1, 28/11.24]), 1, [1, 2, 3]);
        fm.addFeature(SynEM.FeatureLegacy.PrincipalAxis('prod', ...
            [1, 1, 28/11.24]), [2 3], 1);
        fm.addFeature(SynEM.FeatureLegacy.ConvexHull(3),[1, 2, 3],[1 3 2]);

    case 'CnnOnly'
        %fms from last two layers
        
        %check for cnet
        if isempty(varargin) || ~isa(varargin{1}, 'Codat.CNN.cnn')
            error(['Missing input: Specify a Codat.CNN object in ' ...
                'varargin{1}.'])
        end
        cnet = varargin{1};
        
        %varargin{2} specifies which fms of the cnn to use (default last
        %hidden layer and output layer)
        if length(varargin) == 2
            fmL = varargin{2};
        else
            fmL = [cnet.layer-1, cnet.layer];
        end
        
        areaT = 150;
        subvolsSize = [40 80 160];
        quantiles = [0.25 0.5 0.75 0.01 0.99];
        moments = true(4,1);
        voxelSize = [11.24, 11.24, 28];
        fRawNorm = @(x)(single(x) - 122)./22;

        %construct feature map
        fm = FeatureMap(subvolsSize, areaT, ...
            quantiles, moments, 'direction', voxelSize, fRawNorm);

        %add features
        fm.addFeature(Identity());
        fm.addFeature(CNNFeat(cnet, ...
            struct('target_size', [64, 64, 64], 'convAlg', 'fft2', ...
            'fmL', fmL, 'display', 10)));
        fm.addFeature(Volume([], 3), [1, 2, 3], [1 3 2]);
        fm.addFeature(PrincipalAxis('length'), 1 ,[1, 2, 3]);
        fm.addFeature(PrincipalAxis('prod'), [2 3], 1);

    case 'SurfaceOnlyTest'
        areaT = 150;
        subvolsSize = [];
        quantiles = [0.25 0.5 0.75 0 1];
        moments = true(4,1);
        voxelSize = [11.24, 11.24, 28];
        sigma = 12./voxelSize;

        %construct feature map
        fm = FeatureMap(subvolsSize, areaT, ...
            quantiles, moments, 'direction');

        %add features
        fm.addFeature(Identity());
        fm.addFeature(GaussFilter(sigma, 2));
        fm.addFeature(Volume([], 1),1,1);

    otherwise
        error('Feature map %s not defined.', name);
end

%select all features for now
fm.setSelectedFeat();

end
