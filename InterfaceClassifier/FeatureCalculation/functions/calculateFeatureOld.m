function X = calculateFeatureOld( name, parameters, raw, interfaceSurfaceList, subsegmentsList, featureIndices, quantiles, sumStats, border )
%CALCULATEFEATURE Calculates the specified feature.
% OUTPUT X: Array, where the rows correspond to interfaces and the columns
%           correspond to different features.
%
% NOTE see makeFeatureMapOld.
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

fprintf('%s INFO - calculating feature %s\n',datestr(now,'yyyy-mm-dd HH:MM:SS,FFF'),name);
switch name
    
    case 'Identity'
        imfeat = raw;
        X = poolFromFilterResponse(imfeat,interfaceSurfaceList,subsegmentsList,featureIndices{1},quantiles,sumStats,border);
        
    case 'GaussianSmoothed'
        imfeat = GaussianSmoothed(raw,parameters{1},parameters{2});
        X = poolFromFilterResponse(imfeat,interfaceSurfaceList,subsegmentsList,featureIndices{1},quantiles,sumStats,border);
        
    case 'LoG'
        imfeat = LaplacianOfGaussian(raw,parameters{1},parameters{2});
        X = poolFromFilterResponse(imfeat,interfaceSurfaceList,subsegmentsList,featureIndices{1},quantiles,sumStats,border);
        
    case 'DoG'
        imfeat = DifferenceOfGaussians(raw,parameters{1},parameters{2},parameters{3});
        X = poolFromFilterResponse(imfeat,interfaceSurfaceList,subsegmentsList,featureIndices{1},quantiles,sumStats,border);
        
    case 'GaussGrad'
        imfeat = GaussGradientMagnitude(raw,parameters{1},parameters{2});
        X = poolFromFilterResponse(imfeat,interfaceSurfaceList,subsegmentsList,featureIndices{1},quantiles,sumStats,border);
        
    case 'StructureTensor'
        imfeat = StructureTensor(raw,parameters{1},parameters{2},parameters{3},parameters{4});
        X1 = poolFromFilterResponse(imfeat{1},interfaceSurfaceList,subsegmentsList,featureIndices{1},quantiles,sumStats,border);
        X2 = poolFromFilterResponse(imfeat{2},interfaceSurfaceList,subsegmentsList,featureIndices{2},quantiles,sumStats,border);
        X3 = poolFromFilterResponse(imfeat{3},interfaceSurfaceList,subsegmentsList,featureIndices{3},quantiles,sumStats,border);
        X = cat(2,X1,X2,X3);
        
    case 'Hessian'
        imfeat = Hessian(raw,parameters{1},parameters{2});
        X1 = poolFromFilterResponse(imfeat{1},interfaceSurfaceList,subsegmentsList,featureIndices{1},quantiles,sumStats,border);
        X2 = poolFromFilterResponse(imfeat{2},interfaceSurfaceList,subsegmentsList,featureIndices{2},quantiles,sumStats,border);
        X3 = poolFromFilterResponse(imfeat{3},interfaceSurfaceList,subsegmentsList,featureIndices{3},quantiles,sumStats,border);
        X = cat(2,X1,X2,X3);
        
    case 'StdDev'
        nhood = parameters{1};
        imfeat = stdfilt(raw,nhood);
        X = poolFromFilterResponse(imfeat,interfaceSurfaceList,subsegmentsList,featureIndices{1},quantiles,sumStats,border);
        
    case 'Entropy'
        nhood = parameters{1};
        imfeat = single(entropyfilt(uint8(raw),nhood));
        X = poolFromFilterResponse(imfeat,interfaceSurfaceList,subsegmentsList,featureIndices{1},quantiles,sumStats,border);
        
    case 'lPol'
        h = parameters{1};
        raw = raw./255;
        imfeat = imfilter(raw.*raw,h) - imfilter(raw,h).^2;
        X = poolFromFilterResponse(imfeat,interfaceSurfaceList,subsegmentsList,featureIndices{1},quantiles,sumStats,border);
        
    case 'BallAverage'
        imfeat = BallAverage(raw,parameters{1});
        X = poolFromFilterResponse(imfeat,interfaceSurfaceList,subsegmentsList,featureIndices{1},quantiles,sumStats,border);
        
    case 'MaximumFilter'
        nhood = parameters{1};
        sigma = parameters{2};
        rawSmoothed = Codat.Features.gaussianFilter(raw,sigma);
        imfeat = Codat.Features.maximumFilter(rawSmoothed,nhood);
        X = poolFromFilterResponse(imfeat,interfaceSurfaceList,subsegmentsList,featureIndices{1},quantiles,sumStats,border);
        
    case 'MinimumFilter'
        nhood = parameters{1};
        sigma = parameters{2};
        rawSmoothed = Codat.Features.gaussianFilter(raw,sigma);
        imfeat = Codat.Features.minimumFilter(rawSmoothed,nhood);
        X = poolFromFilterResponse(imfeat,interfaceSurfaceList,subsegmentsList,featureIndices{1},quantiles,sumStats,border);
        
    case 'Area'
        indices = featureIndices{1};
        volumes = cat(2,interfaceSurfaceList,subsegmentsList(:,end-1:end));
        volumes(:,~indices) = [];
        if isempty(volumes)
            X = [];
        else
            X = zeros(size(volumes,1),size(volumes,2),'single');
            for k = 1:size(volumes,1)
                for l = 1:size(volumes,2)
                    X(k,l) = length(volumes{k,l});
                end
            end
            if size(X,2) == 1
                X = cat(1,X,X);
            elseif size(X,2) == 2
                X = cat(1,X,X(:,[2,1]));
            else
                X = cat(1,X,X(:,[1,3,2]));
            end
        end    
        
    case 'Diameter'
        surfaceDiameter = zeros(length(interfaceSurfaceList),1,'single');
        if featureIndices{1}
            for k = 1:length(interfaceSurfaceList)
                surfaceDiameter(k) = single(nthroot(6*length(interfaceSurfaceList{k})*(11.24*11.24*28)*pi,3));
            end
            X = cat(1,surfaceDiameter,surfaceDiameter);
        else
            X = [];
        end
        
    case 'PrincipalAxis'
        indices = featureIndices{1};
        if any(indices)
            axisLength = zeros(length(interfaceSurfaceList),3,'single');
            for k = 1:length(interfaceSurfaceList)
                [x,y,z] = ind2sub(size(raw),interfaceSurfaceList{k});
                [~,~,latent] = princomp(double([x,y,28/11.24.*z]),'econ');
                %if the contact surface is a plane set third latent to zero
                if length(latent) < 3
                    latent(end + 1:3) = single(0); 
                end
                axisLength(k,:) = single(latent);
            end
            rows = subsref(1:3,struct('type','()','subs',{{indices}}));
            X = cat(1,axisLength,axisLength(:,[1 3 2]));
            X = X(:,rows);
        else
            X = [];
        end
        
    case 'PAxisProduct'
        majorAxisProd = zeros(size(subsegmentsList,1),1,'single');
        if featureIndices{1}
            for k = 1:size(subsegmentsList,1)
                [x1,y1,z1] = ind2sub(size(raw),subsegmentsList{k,end-1});
                [x2,y2,z2] = ind2sub(size(raw),subsegmentsList{k,end});
                if length(x1) < 2 || length(x2) < 2
                    majorAxisProd(k) = 0;
                else
                    try
                        pca1 = princomp([x1,y1,28/11.24.*z1],'econ');
                    catch
                        C = cov(single([x1,y1,28/11.24.*z1]));
                        pca1 = pcacov(C);
                    end
                    try
                        pca2 = princomp([x2,y2,28/11.24.*z2],'econ');
                    catch
                        C = cov(single([x2,y2,28/11.24.*z2]));
                        pca2 = pcacov(C);
                    end
                    majorAxisProd(k) = pca1(:,1)'*pca2(:,1);
                end
            end
            X = cat(1,majorAxisProd,majorAxisProd);
        else
            X = [];
        end
        
    case 'ConvexHull'
        indices = featureIndices{1};
        volumes = cat(2,interfaceSurfaceList,subsegmentsList(:,end-1:end));
        volumes(:,~indices) = [];
        if isempty(volumes)
            X = [];
        else
            X = zeros(length(interfaceSurfaceList),sum(indices),'single');
            for k = 1:size(volumes,1)
                for l = 1:size(volumes,2)
                    [x,y,z] = ind2sub(size(raw),double(interfaceSurfaceList{k}));
                    try
                        [~,V] = convhull(x,y,z);
                    catch err %if points are colinear or coplanar
                        if length(unique(z)) == 1
                            [~,V] = convhull(x,y);
                        elseif length(unique(x)) == 1
                            [~,V] = convhull(y,z);
                        elseif length(unique(y)) == 1
                            [~,V] = convhull(x,z);
                        else
                            rethrow(err);
                        end
                    end
                    X(k,l) = single(V);
                end
            end
            if size(X,2) == 1
                X = cat(1,X,X);
            elseif size(X,2) == 2
                X = cat(1,X,X(:,[2,1]));
            else
                X = cat(1,X,X(:,[1,3,2]));
            end
        end

    otherwise
        error('Unknown feature specified');
end




end

