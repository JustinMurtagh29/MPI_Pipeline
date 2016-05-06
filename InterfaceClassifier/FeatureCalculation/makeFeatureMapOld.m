function featureMap = makeFeatureMapOld( rinclude,quantiles,summaryStatistics,sumStatNames,areaThreshold, voxelSize)
%MAKEFEATUREMAPOLD Defines features for interface classifier.
% Defines a struct of features containing the information to call the
% calculateFeature function.
% INPUT rinclude: array containing the pooling volume sizes
%                 (Default: [40 80 160]
%       quantiles: array containing the quantiles
%                  (Default: [0.25,0.5,0.75])
%       summaryStatistics: cell array containing the names/abbreviations of
%                  the summaryStatistics
%                  (Default: see sumStatNames)
%       sumStatNames: Names/abbreviations of summary statistics
%                  (Default: {'min','max','mean','var','skew','kur'})
%       areaThreshold: Interface area threshold
%                  (Default: 150)
%       voxelSize: [1x3] array of double containing the voxel size in nm.
%                  (Default: [11.24, 11.24, 28])
% OUTPUT featureMap: struct containing the names and parameters of the
%                    features required for the call of calculateFeatures
%                    and calculateInterfacesFromSeg
%
% NOTE This version can be used together with calculateFeaturesOld. To do
%      so the calculateFeature function in calculateFeatures needs to be
%      changes to calculateFeaturesOld.
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%set default settings
if ~exist('rinclude','var') || isempty(rinclude)
    rinclude = [40 80 160];
end
if ~exist('quantiles','var') || isempty(quantiles)
    quantiles = [0.25,0.5,0.75];
end
if ~exist('summaryStatistics','var') || isempty(summaryStatistics)
    summaryStatistics = cellfun(@func2str,{@(x)single(min(x)), ...
        @(x)single(max(x)),@(x)mean(x),@(x)var(x),@(x)skewness(x), ...
        @(x)kurtosis(x)},'UniformOutput',false);
end
if ~exist('sumStatNames','var') || isempty(sumStatNames)
    sumStatNames = {'min','max','mean','var','skew','kur'};
end
if ~exist('areaThreshold','var') || isempty(areaThreshold)
    areaThreshold = 150;
end
if ~exist('voxelSize','var') || isempty(voxelSize)
    voxelSize = [11.24 11.24 28];
end

%reference sigma and size of filters for pooling
sigma = 12./voxelSize;
filterSiz = ceil(2*sigma);
numQuantiles = length(quantiles);
numSumStat = length(summaryStatistics);

features = { ...
    {'Identity','texture',1,'id'}
    {'StructureTensor','texture',3,'ST1',sigma,filterSiz,sigma,filterSiz}
    {'StructureTensor','texture',3,'ST2',sigma,filterSiz,2.*sigma,2.*filterSiz}
    {'StructureTensor','texture',3,'ST3',2.*sigma,2.*filterSiz,sigma,filterSiz}
    {'StructureTensor','texture',3,'ST4',2.*sigma,2.*filterSiz,2.*sigma,2.*filterSiz}
    {'StructureTensor','texture',3,'ST5',3.*sigma,3.*filterSiz,3.*sigma,3.*filterSiz}
    {'Hessian','texture',3,'H1',sigma,filterSiz};
    {'Hessian','texture',3,'H2',2.*sigma,2.*filterSiz};
    {'Hessian','texture',3,'H3',3.*sigma,3.*filterSiz};
    {'Hessian','texture',3,'H4',4.*sigma,4.*filterSiz};
    {'GaussianSmoothed','texture',1,'Gauss1',sigma,filterSiz};
    {'GaussianSmoothed','texture',1,'Gauss2',2.*sigma,2.*filterSiz};
    {'GaussianSmoothed','texture',1,'Gauss3',3.*sigma,3.*filterSiz};
    {'DoG','texture',1,'DoG1',sigma,filterSiz,1.5};
    {'DoG','texture',1,'DoG2',sigma,filterSiz,2};
    {'DoG','texture',1,'DoG3',2.*sigma,2.*filterSiz,1.5};
    {'DoG','texture',1,'DoG4',2.*sigma,2.*filterSiz,2};
    {'DoG','texture',1,'DoG5',3.*sigma,3.*filterSiz,1.5};
    {'LoG','texture',1,'LoG1',sigma,filterSiz};
    {'LoG','texture',1,'LoG2',2.*sigma,2.*filterSiz};
    {'LoG','texture',1,'LoG3',3.*sigma,3.*filterSiz};
    {'LoG','texture',1,'LoG4',4.*sigma,4.*filterSiz};
    {'GaussGrad','texture',1,'GaussGrad1',sigma,filterSiz};
    {'GaussGrad','texture',1,'GaussGrad2',2.*sigma,2.*filterSiz};
    {'GaussGrad','texture',1,'GaussGrad3',3.*sigma,3.*filterSiz};
    {'GaussGrad','texture',1,'GaussGrad4',4.*sigma,4.*filterSiz};
    {'GaussGrad','texture',1,'GaussGrad5',5.*sigma,5.*filterSiz};
    {'StdDev','texture',1,'StdDev',ones(5,5,5)};
    {'Entropy','texture',1,'Entropy',ones(5,5,5)};
    {'lPol','texture',1,'lPol3',ones(3,3,3)};
    {'lPol','texture',1,'lPol5',ones(5,5,5)};
    {'BallAverage','texture',1,'BallAverage3',3};
    {'BallAverage','texture',1,'BallAverage6',6};
    {'Area','shape',3,'area'};
    {'Diameter','shape',1,'diameter'};
    {'PrincipalAxis','shape',3,'paxis'};
    {'PAxisProduct','shape',1,'paxisprod'};
    {'ConvexHull','shape',3,'convhull'};
    };

names = cellfun(@(x)x{1},features,'UniformOutput',false);
type = cellfun(@(x)x{2},features,'UniformOutput',false);
multiplier = cell2mat(cellfun(@(x)x{3},features,'UniformOutput',false));
id = cellfun(@(x)x{4},features,'UniformOutput',false);
parameters = cellfun(@(x)x(5:end),features,'UniformOutput',false);

featureMap = struct;
featureMap.names = names;
featureMap.type = type;
featureMap.multiplier = multiplier;
featureMap.id = id;
featureMap.parameters = parameters;
featureMap.quantiles = quantiles;
featureMap.summaryStatistics = summaryStatistics;
featureMap.summaryStatisticsNames = sumStatNames;
featureMap.areaThreshold = areaThreshold;
featureMap.voxelSize = voxelSize;

isTexture = cellfun(@(x)strcmp(x,'texture'),type);
featureMap.numFeatures = zeros(length(id),1);
featureMap.numFeatures(isTexture) = multiplier(isTexture).*(1+2*length(rinclude)).*(numQuantiles + numSumStat);
featureMap.numFeatures(~isTexture) = multiplier(~isTexture);
featureMap.numSelectedFeatures = sum(featureMap.numFeatures);

%create feature info cell array
featureInfo = cell(0);
endingStrings = getEndingStrings(rinclude,quantiles,sumStatNames);
for k = 1:length(features)
    for l = 1:multiplier(k)
        switch type{k}
            case 'texture'
                featureInfo = cat(1,featureInfo,strcat(id{k},endingStrings,'_f',num2str(l)));
            case 'shape'
                featureInfo = cat(1,featureInfo,strcat(id{k},'_f',num2str(l)));
            otherwise
                error('Unrecognized feature type');
        end
    end
end
featureMap.featureInfo = featureInfo;
featureMap.rinclude = rinclude;
featureMap = getFeatureIndicesOld( featureMap );

end

function endingStrings = getEndingStrings(rinclude,quantiles,sumStat)
quantileNames = cell(length(quantiles),1);
for k = 1:length(quantiles)
    quantileNames{k} = sprintf('q%d',k);
end
poolingStrings = cat(1,quantileNames,sumStat');
poolingStrings = strcat('_',poolingStrings);
endingStrings = strcat(poolingStrings,'_c_r0');
for k = 1:length(rinclude)
    for l = 1:2
        endingStrings = cat(1,endingStrings,strcat(poolingStrings,{sprintf('_s%d_r%d',l,rinclude(k))}));
    end
end
end
