% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

% See https://gitlab.mpcdf.mpg.de/connectomics/pipeline/blob/2f7e9213809260349effed7379c90dea3b35387c/+connectEM/predictLocalCube.m
connectEmFile = '/gaba/u/mberning/results/edgeClassifier/20170729T131244.mat';

% See https://gitlab.mpcdf.mpg.de/connectomics/pipeline/blob/2f7e9213809260349effed7379c90dea3b35387c/+connectEM/classifierTraining.m
connectEmFeatureMapFile = fullfile(rootDir, 'SynapseClassifier.bkp');

info = Util.runInfo();

synEmFile = fullfile( ...
    info.git_repos{1}.local, '+SynEM', ...
    'data', 'SynEMPaperClassifier.mat');

Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);
segSizes = Seg.Global.getSegToSizeMap(param);

%% Load classifiers and feature maps
classifiers = struct;

% SynEm
curClassifier = load(synEmFile);
curClassifier = curClassifier.classifier;

classifiers.SynEm = struct;
classifiers.SynEm.featureMask = ...
    curClassifier.ens.predictorImportance() > 0;
classifiers.SynEm.featureNames = ...
    curClassifier.options.fm.getSelectedFeatureNames();
clear cur*;

% ConnectEM
curClassifier = load(connectEmFile);
curClassifier = curClassifier.classifier;

curFeatureMap = load(connectEmFeatureMapFile, '-mat');
curFeatureMap = curFeatureMap.fm;

classifiers.ConnectEm = struct;
classifiers.ConnectEm.featureMask = ...
    curClassifier.predictorImportance() > 0;
classifiers.ConnectEm.featureNames = repmat( ...
    curFeatureMap.getSelectedFeatureNames(), 2, 1);

%% Segments
segVolDims = 1 + diff(param.bbox, 1, 2);
segVolDims = segVolDims .* param.raw.voxelSize(:) / 1E3 %#ok
segVol = prod(segVolDims) %#ok

% Number of segments
numSegments = maxSegId %#ok

% Segment volume
vxUm3 = prod(param.raw.voxelSize / 1E3);
meanSegSizesUm3 = mean(vxUm3 .* segSizes) %#ok
stdSegSizesUm3 = std(vxUm3 .* segSizes) %#ok

%% SynEM / connectEM features
classifierNames = fieldnames(classifiers);
for curClassifierIdx = 1:numel(classifierNames)
    curClassifierName = classifierNames{curClassifierIdx};
    curStruct = classifiers.(curClassifierName);
    
    curFeatureNames = curStruct.featureNames;
    curFeatureNames = curFeatureNames(curStruct.featureMask);

    first = @(vals) vals(1);
    curFilterNames = cellfun( ....
        @(n) first(strsplit(n, '_')), ...
        curFeatureNames);

    curFilterT = table;
    [curFilterT.names, ~, curUniRows] = unique(curFilterNames);
    curFilterT.isShape = ismember( ...
        curFilterT.names, {'Volume', 'PrAx', 'ConvHull'});
    curFilterT.numFeatures = accumarray(curUniRows, 1);
    
    disp(curClassifierName)
    disp(curFilterT)
    
    curNumShapeFeatures = sum( ...
        curFilterT.numFeatures(curFilterT.isShape)) %#ok
    curNumTextureFeatures = sum( ...
        curFilterT.numFeatures(~curFilterT.isShape)) %#ok
end
