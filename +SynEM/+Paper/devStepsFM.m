function fm = devStepsFM( step )
%DEVSTEPSFM Feature maps for dev steps.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

import SynEM.Feature.*
import SynEM.FeatureMap


switch step
    case 1
        %no moments, no shape features (except volume), no dog, stdDev,
        %entropy lPol, ballAverage, no direction, only 160 nm subvol
        areaT = 150;
        subvolsSize = [160];
        quantiles = [0.25 0.5 0.75 0 1];
        voxelSize = [11.24, 11.24, 28];
        moments = false(4,1);
        sigma = 12./voxelSize;
        fS = ceil(2.*sigma);
        fRawNorm = '@(x)single(x)';

        %construct feature map
        fm = FeatureMap(subvolsSize, areaT, ...
            quantiles, moments, 'single', [1, 1, 1], fRawNorm);


        %add features
        fm.addFeature(Identity());
        fm.addFeature(EVsStructureTensor(sigma, fS, sigma, fS));
        fm.addFeature(EVsStructureTensor(sigma, fS, 2.*sigma, 2.*fS));
        fm.addFeature(EVsStructureTensor(2.*sigma, 2.*fS, sigma, fS));
        fm.addFeature(EVsStructureTensor(2.*sigma, 2.*fS, 2.*sigma, 2.*fS));
        fm.addFeature(EVsStructureTensor(3.*sigma, 3.*fS, 3.*sigma, 3.*fS));
        fm.addFeature(EVsHessian(sigma, fS));
        fm.addFeature(EVsHessian(2.*sigma, 2.*fS));
        fm.addFeature(EVsHessian(3.*sigma, 3.*fS));
        fm.addFeature(EVsHessian(4.*sigma, 4.*fS));
        fm.addFeature(GaussFilter(sigma, fS));
        fm.addFeature(GaussFilter(2.*sigma, 2.*fS));
        fm.addFeature(GaussFilter(3.*sigma, 3.*fS));
        fm.addFeature(LoG(sigma, fS));
        fm.addFeature(LoG(2.*sigma, 2.*fS));
        fm.addFeature(LoG(3.*sigma, 3.*fS));
        fm.addFeature(LoG(4.*sigma, 4.*fS));
        fm.addFeature(GaussGradMagnitude(sigma, fS));
        fm.addFeature(GaussGradMagnitude(2.*sigma, 2.*fS));
        fm.addFeature(GaussGradMagnitude(3.*sigma, 3.*fS));
        fm.addFeature(GaussGradMagnitude(4.*sigma, 4.*fS));
        fm.addFeature(GaussGradMagnitude(5.*sigma, 5.*fS));
        fm.addFeature(Volume([], 3),[1, 2, 3],[1 3 2]);

    case 2
        %no moments, no direction, only 160 nm subvol
        areaT = 150;
        subvolsSize = [160];
        quantiles = [0.25 0.5 0.75 0 1];
        voxelSize = [11.24, 11.24, 28];
        moments = false(4,1);
        sigma = 12./voxelSize;
        fS = ceil(2.*sigma);
        fRawNorm = '@(x)single(x)';

        %construct feature map
        fm = FeatureMap(subvolsSize, areaT, ...
            quantiles, moments, 'single', [1, 1, 1], fRawNorm);

        %add features
        fm.addFeature(Identity());
        fm.addFeature(EVsStructureTensor(sigma, fS, sigma, fS));
        fm.addFeature(EVsStructureTensor(sigma, fS, 2.*sigma, 2.*fS));
        fm.addFeature(EVsStructureTensor(2.*sigma, 2.*fS, sigma, fS));
        fm.addFeature(EVsStructureTensor(2.*sigma, 2.*fS, 2.*sigma, 2.*fS));
        fm.addFeature(EVsStructureTensor(3.*sigma, 3.*fS, 3.*sigma, 3.*fS));
        fm.addFeature(EVsHessian(sigma, fS));
        fm.addFeature(EVsHessian(2.*sigma, 2.*fS));
        fm.addFeature(EVsHessian(3.*sigma, 3.*fS));
        fm.addFeature(EVsHessian(4.*sigma, 4.*fS));
        fm.addFeature(GaussFilter(sigma, fS));
        fm.addFeature(GaussFilter(2.*sigma, 2.*fS));
        fm.addFeature(GaussFilter(3.*sigma, 3.*fS));
        fm.addFeature(DoG(sigma, fS, 1.5));
        fm.addFeature(DoG(sigma, fS, 2));
        fm.addFeature(DoG(2.*sigma, 2.*fS, 1.5));
        fm.addFeature(DoG(2.*sigma, 2.*fS, 2));
        fm.addFeature(DoG(3.*sigma, 3.*fS, 1.5));
        fm.addFeature(LoG(sigma, fS));
        fm.addFeature(LoG(2.*sigma, 2.*fS));
        fm.addFeature(LoG(3.*sigma, 3.*fS));
        fm.addFeature(LoG(4.*sigma, 4.*fS));
        fm.addFeature(GaussGradMagnitude(sigma, fS));
        fm.addFeature(GaussGradMagnitude(2.*sigma, 2.*fS));
        fm.addFeature(GaussGradMagnitude(3.*sigma, 3.*fS));
        fm.addFeature(GaussGradMagnitude(4.*sigma, 4.*fS));
        fm.addFeature(GaussGradMagnitude(5.*sigma, 5.*fS));
        fm.addFeature(StdFilter([5, 5, 5]));
        fm.addFeature(EntropyFilter([5, 5, 5]));
        fm.addFeature(IntVar([3, 3, 3]));
        fm.addFeature(IntVar([5, 5, 5]));
        fm.addFeature(AverageFilter('ball',3,[1, 1, 2]));
        fm.addFeature(AverageFilter('ball',6,[1, 1, 2]));
        fm.addFeature(Volume([], 3),[1, 2, 3],[1 3 2]);
        fm.addFeature(Volume( ...
            '@(x)single(nthroot(6*x*(11.24*11.24*28)*pi,3))',1), ...
            1, 1); %diameter from old fm
        fm.addFeature(PrincipalAxis('length',[1, 1, 28/11.24]),1,[1, 2, 3]);
        fm.addFeature(PrincipalAxis('prod',[1, 1, 28/11.24]),[2 3],1);
        fm.addFeature(ConvexHull(3),[1, 2, 3],[1 3 2]);

    case 3
        %no moments, no direction
        areaT = 150;
        subvolsSize = [40 80 160];
        quantiles = [0.25 0.5 0.75 0 1];
        voxelSize = [11.24, 11.24, 28];
        moments = false(4,1);
        sigma = 12./voxelSize;
        fS = ceil(2.*sigma);
        fRawNorm = '@(x)single(x)';

        %construct feature map
        fm = FeatureMap(subvolsSize, areaT, ...
            quantiles, moments, 'single', [1, 1, 1], fRawNorm);

        %add features
        fm.addFeature(Identity());
        fm.addFeature(EVsStructureTensor(sigma, fS, sigma, fS));
        fm.addFeature(EVsStructureTensor(sigma, fS, 2.*sigma, 2.*fS));
        fm.addFeature(EVsStructureTensor(2.*sigma, 2.*fS, sigma, fS));
        fm.addFeature(EVsStructureTensor(2.*sigma, 2.*fS, 2.*sigma, 2.*fS));
        fm.addFeature(EVsStructureTensor(3.*sigma, 3.*fS, 3.*sigma, 3.*fS));
        fm.addFeature(EVsHessian(sigma, fS));
        fm.addFeature(EVsHessian(2.*sigma, 2.*fS));
        fm.addFeature(EVsHessian(3.*sigma, 3.*fS));
        fm.addFeature(EVsHessian(4.*sigma, 4.*fS));
        fm.addFeature(GaussFilter(sigma, fS));
        fm.addFeature(GaussFilter(2.*sigma, 2.*fS));
        fm.addFeature(GaussFilter(3.*sigma, 3.*fS));
        fm.addFeature(DoG(sigma, fS, 1.5));
        fm.addFeature(DoG(sigma, fS, 2));
        fm.addFeature(DoG(2.*sigma, 2.*fS, 1.5));
        fm.addFeature(DoG(2.*sigma, 2.*fS, 2));
        fm.addFeature(DoG(3.*sigma, 3.*fS, 1.5));
        fm.addFeature(LoG(sigma, fS));
        fm.addFeature(LoG(2.*sigma, 2.*fS));
        fm.addFeature(LoG(3.*sigma, 3.*fS));
        fm.addFeature(LoG(4.*sigma, 4.*fS));
        fm.addFeature(GaussGradMagnitude(sigma, fS));
        fm.addFeature(GaussGradMagnitude(2.*sigma, 2.*fS));
        fm.addFeature(GaussGradMagnitude(3.*sigma, 3.*fS));
        fm.addFeature(GaussGradMagnitude(4.*sigma, 4.*fS));
        fm.addFeature(GaussGradMagnitude(5.*sigma, 5.*fS));
        fm.addFeature(StdFilter([5, 5, 5]));
        fm.addFeature(EntropyFilter([5, 5, 5]));
        fm.addFeature(IntVar([3, 3, 3]));
        fm.addFeature(IntVar([5, 5, 5]));
        fm.addFeature(AverageFilter('ball',3,[1, 1, 2]));
        fm.addFeature(AverageFilter('ball',6,[1, 1, 2]));
        fm.addFeature(Volume([], 3),[1, 2, 3],[1 3 2]);
        fm.addFeature(Volume( ...
            '@(x)single(nthroot(6*x*(11.24*11.24*28)*pi,3))',1), ...
            1, 1); %diameter from old fm
        fm.addFeature(PrincipalAxis('length',[1, 1, 28/11.24]),1,[1, 2, 3]);
        fm.addFeature(PrincipalAxis('prod',[1, 1, 28/11.24]),[2 3],1);
        fm.addFeature(ConvexHull(3),[1, 2, 3],[1 3 2]);

    case 4
        %no direction
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
            quantiles, moments, 'single', [1, 1, 1], fRawNorm);

        %add features
        fm.addFeature(Identity());
        fm.addFeature(EVsStructureTensor(sigma, fS, sigma, fS));
        fm.addFeature(EVsStructureTensor(sigma, fS, 2.*sigma, 2.*fS));
        fm.addFeature(EVsStructureTensor(2.*sigma, 2.*fS, sigma, fS));
        fm.addFeature(EVsStructureTensor(2.*sigma, 2.*fS, 2.*sigma, 2.*fS));
        fm.addFeature(EVsStructureTensor(3.*sigma, 3.*fS, 3.*sigma, 3.*fS));
        fm.addFeature(EVsHessian(sigma, fS));
        fm.addFeature(EVsHessian(2.*sigma, 2.*fS));
        fm.addFeature(EVsHessian(3.*sigma, 3.*fS));
        fm.addFeature(EVsHessian(4.*sigma, 4.*fS));
        fm.addFeature(GaussFilter(sigma, fS));
        fm.addFeature(GaussFilter(2.*sigma, 2.*fS));
        fm.addFeature(GaussFilter(3.*sigma, 3.*fS));
        fm.addFeature(DoG(sigma, fS, 1.5));
        fm.addFeature(DoG(sigma, fS, 2));
        fm.addFeature(DoG(2.*sigma, 2.*fS, 1.5));
        fm.addFeature(DoG(2.*sigma, 2.*fS, 2));
        fm.addFeature(DoG(3.*sigma, 3.*fS, 1.5));
        fm.addFeature(LoG(sigma, fS));
        fm.addFeature(LoG(2.*sigma, 2.*fS));
        fm.addFeature(LoG(3.*sigma, 3.*fS));
        fm.addFeature(LoG(4.*sigma, 4.*fS));
        fm.addFeature(GaussGradMagnitude(sigma, fS));
        fm.addFeature(GaussGradMagnitude(2.*sigma, 2.*fS));
        fm.addFeature(GaussGradMagnitude(3.*sigma, 3.*fS));
        fm.addFeature(GaussGradMagnitude(4.*sigma, 4.*fS));
        fm.addFeature(GaussGradMagnitude(5.*sigma, 5.*fS));
        fm.addFeature(StdFilter([5, 5, 5]));
        fm.addFeature(EntropyFilter([5, 5, 5]));
        fm.addFeature(IntVar([3, 3, 3]));
        fm.addFeature(IntVar([5, 5, 5]));
        fm.addFeature(AverageFilter('ball',3,[1, 1, 2]));
        fm.addFeature(AverageFilter('ball',6,[1, 1, 2]));
        fm.addFeature(Volume([], 3),[1, 2, 3],[1 3 2]);
        fm.addFeature(Volume( ...
            '@(x)single(nthroot(6*x*(11.24*11.24*28)*pi,3))',1), ...
            1, 1); %diameter from old fm
        fm.addFeature(PrincipalAxis('length',[1, 1, 28/11.24]),1,[1, 2, 3]);
        fm.addFeature(PrincipalAxis('prod',[1, 1, 28/11.24]),[2 3],1);
        fm.addFeature(ConvexHull(3),[1, 2, 3],[1 3 2]);

    case {5, 6}
        %final feature map
        % step 5 with ada boost
        % step 6 logit boost
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
        fm.addFeature(Identity());
        fm.addFeature(EVsStructureTensor(sigma, fS, sigma, fS));
        fm.addFeature(EVsStructureTensor(sigma, fS, 2.*sigma, 2.*fS));
        fm.addFeature(EVsStructureTensor(2.*sigma, 2.*fS, sigma, fS));
        fm.addFeature(EVsStructureTensor(2.*sigma, 2.*fS, 2.*sigma, 2.*fS));
        fm.addFeature(EVsStructureTensor(3.*sigma, 3.*fS, 3.*sigma, 3.*fS));
        fm.addFeature(EVsHessian(sigma, fS));
        fm.addFeature(EVsHessian(2.*sigma, 2.*fS));
        fm.addFeature(EVsHessian(3.*sigma, 3.*fS));
        fm.addFeature(EVsHessian(4.*sigma, 4.*fS));
        fm.addFeature(GaussFilter(sigma, fS));
        fm.addFeature(GaussFilter(2.*sigma, 2.*fS));
        fm.addFeature(GaussFilter(3.*sigma, 3.*fS));
        fm.addFeature(DoG(sigma, fS, 1.5));
        fm.addFeature(DoG(sigma, fS, 2));
        fm.addFeature(DoG(2.*sigma, 2.*fS, 1.5));
        fm.addFeature(DoG(2.*sigma, 2.*fS, 2));
        fm.addFeature(DoG(3.*sigma, 3.*fS, 1.5));
        fm.addFeature(LoG(sigma, fS));
        fm.addFeature(LoG(2.*sigma, 2.*fS));
        fm.addFeature(LoG(3.*sigma, 3.*fS));
        fm.addFeature(LoG(4.*sigma, 4.*fS));
        fm.addFeature(GaussGradMagnitude(sigma, fS));
        fm.addFeature(GaussGradMagnitude(2.*sigma, 2.*fS));
        fm.addFeature(GaussGradMagnitude(3.*sigma, 3.*fS));
        fm.addFeature(GaussGradMagnitude(4.*sigma, 4.*fS));
        fm.addFeature(GaussGradMagnitude(5.*sigma, 5.*fS));
        fm.addFeature(StdFilter([5, 5, 5]));
        fm.addFeature(EntropyFilter([5, 5, 5]));
        fm.addFeature(IntVar([3, 3, 3]));
        fm.addFeature(IntVar([5, 5, 5]));
        fm.addFeature(AverageFilter('ball',3,[1, 1, 2]));
        fm.addFeature(AverageFilter('ball',6,[1, 1, 2]));
        fm.addFeature(Volume([], 3),[1, 2, 3],[1 3 2]);
        fm.addFeature(Volume( ...
            '@(x)single(nthroot(6*x*(11.24*11.24*28)*pi,3))',1), ...
            1, 1); %diameter from old fm
        fm.addFeature(PrincipalAxis('length',[1, 1, 28/11.24]),1,[1, 2, 3]);
        fm.addFeature(PrincipalAxis('prod',[1, 1, 28/11.24]),[2 3],1);
        fm.addFeature(ConvexHull(3),[1, 2, 3],[1 3 2]);

    case 7
        %basically final feature map but with single (data augmentation
        %test)
        % step 7 augment label
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
        fm.addFeature(Identity());
        fm.addFeature(EVsStructureTensor(sigma, fS, sigma, fS));
        fm.addFeature(EVsStructureTensor(sigma, fS, 2.*sigma, 2.*fS));
        fm.addFeature(EVsStructureTensor(2.*sigma, 2.*fS, sigma, fS));
        fm.addFeature(EVsStructureTensor(2.*sigma, 2.*fS, 2.*sigma, 2.*fS));
        fm.addFeature(EVsStructureTensor(3.*sigma, 3.*fS, 3.*sigma, 3.*fS));
        fm.addFeature(EVsHessian(sigma, fS));
        fm.addFeature(EVsHessian(2.*sigma, 2.*fS));
        fm.addFeature(EVsHessian(3.*sigma, 3.*fS));
        fm.addFeature(EVsHessian(4.*sigma, 4.*fS));
        fm.addFeature(GaussFilter(sigma, fS));
        fm.addFeature(GaussFilter(2.*sigma, 2.*fS));
        fm.addFeature(GaussFilter(3.*sigma, 3.*fS));
        fm.addFeature(DoG(sigma, fS, 1.5));
        fm.addFeature(DoG(sigma, fS, 2));
        fm.addFeature(DoG(2.*sigma, 2.*fS, 1.5));
        fm.addFeature(DoG(2.*sigma, 2.*fS, 2));
        fm.addFeature(DoG(3.*sigma, 3.*fS, 1.5));
        fm.addFeature(LoG(sigma, fS));
        fm.addFeature(LoG(2.*sigma, 2.*fS));
        fm.addFeature(LoG(3.*sigma, 3.*fS));
        fm.addFeature(LoG(4.*sigma, 4.*fS));
        fm.addFeature(GaussGradMagnitude(sigma, fS));
        fm.addFeature(GaussGradMagnitude(2.*sigma, 2.*fS));
        fm.addFeature(GaussGradMagnitude(3.*sigma, 3.*fS));
        fm.addFeature(GaussGradMagnitude(4.*sigma, 4.*fS));
        fm.addFeature(GaussGradMagnitude(5.*sigma, 5.*fS));
        fm.addFeature(StdFilter([5, 5, 5]));
        fm.addFeature(EntropyFilter([5, 5, 5]));
        fm.addFeature(IntVar([3, 3, 3]));
        fm.addFeature(IntVar([5, 5, 5]));
        fm.addFeature(AverageFilter('ball',3,[1, 1, 2]));
        fm.addFeature(AverageFilter('ball',6,[1, 1, 2]));
        fm.addFeature(Volume([], 3),[1, 2, 3],[1 3 2]);
        fm.addFeature(Volume( ...
            '@(x)single(nthroot(6*x*(11.24*11.24*28)*pi,3))',1), ...
            1, 1); %diameter from old fm
        fm.addFeature(PrincipalAxis('length',[1, 1, 28/11.24]),1,[1, 2, 3]);
        fm.addFeature(PrincipalAxis('prod',[1, 1, 28/11.24]),[2 3],1);
        fm.addFeature(ConvexHull(3),[1, 2, 3],[1 3 2]);

    case 8
        %basically final feature map but with single (data augmentation
        %test)
        % step 8 using single label
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
            quantiles, moments, 'single', [1, 1, 1], fRawNorm);

        %add features
        fm.addFeature(Identity());
        fm.addFeature(EVsStructureTensor(sigma, fS, sigma, fS));
        fm.addFeature(EVsStructureTensor(sigma, fS, 2.*sigma, 2.*fS));
        fm.addFeature(EVsStructureTensor(2.*sigma, 2.*fS, sigma, fS));
        fm.addFeature(EVsStructureTensor(2.*sigma, 2.*fS, 2.*sigma, 2.*fS));
        fm.addFeature(EVsStructureTensor(3.*sigma, 3.*fS, 3.*sigma, 3.*fS));
        fm.addFeature(EVsHessian(sigma, fS));
        fm.addFeature(EVsHessian(2.*sigma, 2.*fS));
        fm.addFeature(EVsHessian(3.*sigma, 3.*fS));
        fm.addFeature(EVsHessian(4.*sigma, 4.*fS));
        fm.addFeature(GaussFilter(sigma, fS));
        fm.addFeature(GaussFilter(2.*sigma, 2.*fS));
        fm.addFeature(GaussFilter(3.*sigma, 3.*fS));
        fm.addFeature(DoG(sigma, fS, 1.5));
        fm.addFeature(DoG(sigma, fS, 2));
        fm.addFeature(DoG(2.*sigma, 2.*fS, 1.5));
        fm.addFeature(DoG(2.*sigma, 2.*fS, 2));
        fm.addFeature(DoG(3.*sigma, 3.*fS, 1.5));
        fm.addFeature(LoG(sigma, fS));
        fm.addFeature(LoG(2.*sigma, 2.*fS));
        fm.addFeature(LoG(3.*sigma, 3.*fS));
        fm.addFeature(LoG(4.*sigma, 4.*fS));
        fm.addFeature(GaussGradMagnitude(sigma, fS));
        fm.addFeature(GaussGradMagnitude(2.*sigma, 2.*fS));
        fm.addFeature(GaussGradMagnitude(3.*sigma, 3.*fS));
        fm.addFeature(GaussGradMagnitude(4.*sigma, 4.*fS));
        fm.addFeature(GaussGradMagnitude(5.*sigma, 5.*fS));
        fm.addFeature(StdFilter([5, 5, 5]));
        fm.addFeature(EntropyFilter([5, 5, 5]));
        fm.addFeature(IntVar([3, 3, 3]));
        fm.addFeature(IntVar([5, 5, 5]));
        fm.addFeature(AverageFilter('ball',3,[1, 1, 2]));
        fm.addFeature(AverageFilter('ball',6,[1, 1, 2]));
        fm.addFeature(Volume([], 3),[1, 2, 3],[1 3 2]);
        fm.addFeature(Volume( ...
            '@(x)single(nthroot(6*x*(11.24*11.24*28)*pi,3))',1), ...
            1, 1); %diameter from old fm
        fm.addFeature(PrincipalAxis('length',[1, 1, 28/11.24]),1,[1, 2, 3]);
        fm.addFeature(PrincipalAxis('prod',[1, 1, 28/11.24]),[2 3],1);
        fm.addFeature(ConvexHull(3),[1, 2, 3],[1 3 2]);

end
