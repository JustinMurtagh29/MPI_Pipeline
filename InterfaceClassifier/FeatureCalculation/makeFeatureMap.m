function featureMap = makeFeatureMap( rinclude,quantiles,summaryStatistics,sumStatNames,areaThreshold,voxelSize)
%MAKEFEATUREMAP Defines features for interface classifier.
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
numQuantiles = length(quantiles);
numSumStat = length(summaryStatistics);

features = { ...
    {'Identity','texture',1,'id'}
%     {'MMV_load','texture',3,'MMV','/path/to/cube'};
    {'StructureTensor','texture',3,'ST1',sigma,sigma}
    {'StructureTensor','texture',3,'ST2',sigma,2.*sigma,2,2}
    {'StructureTensor','texture',3,'ST3',sigma,3.*sigma,2,2}
    {'Hessian','texture',3,'H1',sigma,2};
    {'Hessian','texture',3,'H2',2.*sigma,2};
    {'Hessian','texture',3,'H3',3.*sigma,2};
    {'GaussianSmoothed','texture',1,'Gauss1',sigma};
    {'GaussianSmoothed','texture',1,'Gauss2',2.*sigma,2};
    {'GaussianSmoothed','texture',1,'Gauss3',3.*sigma,2};
    {'DoG','texture',1,'DoG1',sigma,1.5};
    {'DoG','texture',1,'DoG2',sigma,3};
    {'DoG','texture',1,'DoG3',2.*sigma,1.5};
    {'DoG','texture',1,'DoG4',2.*sigma,3};
    {'DoG','texture',1,'DoG5',3.*sigma,1.5,2};
    {'DoG','texture',1,'DoG6',3.*sigma,3,2};
    {'LoG','texture',1,'LoG1',sigma};
    {'LoG','texture',1,'LoG2',2.*sigma};
    {'LoG','texture',1,'LoG3',3.*sigma,2};
    {'LoG','texture',1,'LoG4',4.*sigma,2};
    {'GaussGrad','texture',1,'GaussGrad1',sigma};
    {'GaussGrad','texture',1,'GaussGrad2',2.*sigma,2};
    {'GaussGrad','texture',1,'GaussGrad3',3.*sigma,2};
    {'StdDev','texture',1,'StdDev',[5, 5, 5]};
    {'Entropy','texture',1,'Entropy',[5, 5, 5], sigma};
    {'lPol','texture',1,'lPol3',[3, 3, 3]};
    {'lPol','texture',1,'lPol5',[5, 5, 5]};
    {'SphereAverage','texture',1,'BallAverage3','ball',3};
    {'SphereAverage','texture',1,'BallAverage6','ball',6};
    {'MaximumFilter','texture',1,'maxFilt',3, sigma}
    {'MaximumFilter','texture',1,'maxFilt',5, sigma}
    {'MinimumFilter','texture',1,'minFilt',3, sigma}
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
featureMap = getFeatureIndices( featureMap );

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
