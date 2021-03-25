% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
% Modified by
%   Sahil Loomba <sahil.loomba@brain.mpg.de>

%% Configuration

% Split-biased agglomeration parameters from
% =>20191204T111414aggloGA.mat
% Optimized using
% connectEM.optimizeAgglomeration

% Trying to use this now for glia segmentation.
% NOTE: that segmentAggloPredictions.mat was updated to have 4 classes in June2020 but GA optimization was still done in 2019 with older predictions

%{
param = p;
% m = load(fullfile(param.saveFolder,'aggloGA','20191223T000223optimParams.mat')); % for longer dendrites
m = load(fullfile(param.saveFolder,'aggloGA','20191224T235355optimParams.mat')); % for same axon/dend agglo lengths

aggloParam = struct;
aggloParam.glia.minEdgeProb = m.curParams(1);
aggloParam.glia.minAxonProb = 0;
aggloParam.glia.maxDendProb = 0.5;
aggloParam.glia.maxAstroProb = 1;
aggloParam.glia.minAggloLength = m.curParams(5);
aggloParam.glia.maxVesselProb = 0.5;
aggloParam.glia.maxNucleusProb = 0.5;

datasetNameAppend = ['_ga' '_20191224T235355optimParams'];
timeStamp = [datestr(clock,30) datasetNameAppend];
outDir = fullfile(p.saveFolder, 'aggloMat', [timeStamp '_agglomeration']);

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

parameters.scale.x = num2str(p.raw.voxelSize(1));
parameters.scale.y = num2str(p.raw.voxelSize(2));
parameters.scale.z = num2str(p.raw.voxelSize(3));
parameters.offset.x = '0';
parameters.offset.y = '0';
parameters.offset.z = '0';

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
Util.log('Loading data...')
rootDir = param.saveFolder;
maxSegId = Seg.Global.getMaxSegId(param);

segmentMeta = load([param.saveFolder 'segmentMeta.mat'], 'voxelCount', 'point', 'maxSegId');

segMeta = struct;
segMeta.id = reshape(1:maxSegId, [], 1);
segMeta.pos = Seg.Global.getSegToPointMap(param);

types = load(fullfile(rootDir, 'segmentAggloPredictions.mat'));
heuristics = load(fullfile(rootDir, 'heuristicResult.mat'));

% HACK(amotta): The conversion of classifier outputs to probabilities is
% not perfect. It's thus possible that the probabilities add up to more
% than 100 %. Let's prevent this, so that we can guarantee that dendritic
% and axonal segments will be strictly non-overlapping.
curMask = sum(types.probs, 2) > 1;
types.probs(curMask, :) = ...
    types.probs(curMask, :) ...
 ./ sum(types.probs(curMask, :), 2);

segMeta.probs = nan(maxSegId, numel(types.classes) + 2);
segMeta.probs(types.segId, 1:4) = types.probs;
segMeta.probs(:, 5) = zeros(maxSegId,1); % all non vessel seg have score 0, not NaN
segMeta.probs(heuristics.segIds, 4) = heuristics.vesselScore;
segMeta.probs(:, 6) = zeros(maxSegId,1); % all non nuclei seg have score 0, not NaN
segMeta.probs(heuristics.segIds, 5) = heuristics.nucleiScore;


segMeta.classes = categorical([ ...
    cellstr(types.classes), {'vessel', 'nucleus'}]);
clear types heuristics;

curBorderSize = fullfile(rootDir, 'globalBorder.mat');
curBorderSize = Util.load(curBorderSize, 'borderSize');

% NOTE(amotta): The log transformation counteracts the long tail. This way,
% the border sizes are at least roughly normally distributed.
curBorderSize = log10(double(curBorderSize));
maxBorderSize = max(curBorderSize);

segGraph = load(fullfile(rootDir, 'graph.mat'));
curMask = not(isnan(segGraph.borderIdx));

% NOTE(amotta): Correspondence edges don't have border sizes readily
% available. Let's treat them as infinitely large borders.
segGraph.size = inf(size(segGraph.edges, 1), 1);
segGraph.size(curMask) = curBorderSize(segGraph.borderIdx(curMask));

segGraph = rmfield(segGraph, setdiff( ...
    fieldnames(segGraph), {'edges', 'prob', 'size'}));

%% Reduce multi-graph
Util.log('Doing graph cut...')
clear cur*;

assert(issorted(segGraph.edges, 2));
assert(issortedrows(segGraph.edges));

[~, curSortIds] = sort( ...
    segGraph.prob, 'descend', ...
    'MissingPlacement', 'last');
segGraph.edges = segGraph.edges(curSortIds, :);
segGraph.prob = segGraph.prob(curSortIds);
segGraph.size = segGraph.size(curSortIds);

[~, curKeepIds] = unique(segGraph.edges, 'rows');
segGraph.edges = segGraph.edges(curKeepIds, :);
segGraph.prob = segGraph.prob(curKeepIds);
segGraph.size = segGraph.size(curKeepIds);

%% Build agglomerates
Util.log('Building agglos...')
clear cur*;

[glias, gliaAgglos] = ...
    ConnectEM.buildAgglomerates( ...
        segMeta, segGraph, aggloParam.glia);
%}
curMask = ConnectEM.checkMinLength( ...
    param, aggloParam.glia.minAggloLength, glias);
glias = glias(curMask);
gliaAgglos = gliaAgglos(curMask);

%% Save results
Util.log('Saving results...')
clear cur*;

curOut = struct;
curOut.info = info;
curOut.axons = axons;
curOut.dendrites = dendrites;
curOut.axonAgglos = axonAgglos;
curOut.dendriteAgglos = dendriteAgglos;

curOutFile = sprintf('glia-%s_results.mat', datestr(now, 30));
curOutFile = fullfile(outDir, curOutFile);

Util.saveStruct(curOutFile, curOut);

% Sort agglomerates by voxel size
partition = dendrites;
partitionSize = arrayfun(@(x)sum(segmentMeta.voxelCount(x.segIds)), partition);
[partitionSize, idx] = sort(partitionSize, 'descend');
partition = partition(idx);
dendritesSorted = partition;
dendritesSortedAgglos = arrayfun(@(x) x.segIds,dendritesSorted,'uni',0);

partition = axons;
partitionSize = arrayfun(@(x)sum(segmentMeta.voxelCount(x.segIds)), partition);
[partitionSize, idx] = sort(partitionSize, 'descend');
partition = partition(idx);
axonsSorted = partition;
axonsSortedAgglos = arrayfun(@(x) x.segIds,axonsSorted,'uni',0);

%{
%% write skeletons
Util.log('Export skeletons:')
outputFolderSub = fullfile(outDir,'dendrites');
mkdir(outputFolderSub)
agglosOut = dendritesSortedAgglos(1:100);
display('Writing skeletons for debugging the process:');
parameters.experiment.name = ['Mk1_F6_JS_SubI_v1_mrnet_wsmrnet' '_dend' datasetNameAppend];
Superagglos.skeletonFromAgglo(segGraph.edges, segmentMeta, ...
    agglosOut, 'dendrites', outputFolderSub, parameters);

outputFolderSub = fullfile(outDir,'axons');
mkdir(outputFolderSub)
agglosOut = axonsSortedAgglos(1:100);
display('Writing skeletons for debugging the process:');
parameters.experiment.name = ['Mk1_F6_JS_SubI_v1_mrnet_wsmrnet' '_axon' datasetNameAppend];
Superagglos.skeletonFromAgglo(segGraph.edges, segmentMeta, ...
    agglosOut, 'axons', outputFolderSub, parameters);
%}
maxSegId = segmentMeta.maxSegId; 
%% write segmentation
Util.log('Write new segmentation based on agglos')
segOut = struct;
segOut.root = fullfile(p.saveFolder, 'aggloState', ...
    [timeStamp '_segForGliaAgglos'], '1');
mkdir(segOut.root)
dataType = 'uint32';
wkwInit('new',segOut.root,32, 32, dataType, 1);
segOut.backend = 'wkwrap';
agglosSorted = dendritesSortedAgglos;
mapping = connectEM.createLookup(segmentMeta, agglosSorted);
Seg.Global.applyMappingToSegmentation(p, mapping, segOut);
thisBBox = [1, 1, 1; (ceil(p.bbox(:, 2) ./ 1024) .* 1024)']';
createResolutionPyramid(segOut, thisBBox, [], true);