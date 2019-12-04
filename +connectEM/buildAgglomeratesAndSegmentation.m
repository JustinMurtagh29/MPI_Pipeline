% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
% Modified by
%   Sahil Loomba <sahil.loomba@brain.mpg.de>

%% Configuration

% Split-biased agglomeration parameters from
% =>20191204T111414aggloGA.mat
% Optimized using
% connectEM.optimizeAgglomeration

param = p;
m = load(fullfile(param.saveFolder,'20191204T111414aggloGA.mat'));

aggloParam = struct;
aggloParam.axon.minEdgeProb = m.curParams(1);
aggloParam.axon.minAxonProb = m.curParams(2);
aggloParam.axon.maxDendProb = m.curParams(3);
aggloParam.axon.maxAstroProb = m.curParams(4);
aggloParam.axon.minAggloLength = m.curParams(5);
aggloParam.axon.maxVesselProb = 0.5;
aggloParam.axon.maxNucleusProb = 0.5;

aggloParam.dend.minEdgeProb = m.curParams(6);
aggloParam.dend.minBorderSize = m.curParams(7);
aggloParam.dend.maxAxonProb = m.curParams(8);
aggloParam.dend.minDendProb = m.curParams(9);
aggloParam.dend.maxAstroProb = m.curParams(10);
aggloParam.dend.minAggloLength = m.curParams(11);
aggloParam.dend.maxVesselProb = 0.5;
aggloParam.dend.maxNucleusProb = 0.5;

datasetNameAppend = '_ga';
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
segMeta.probs(types.segId, 1:3) = types.probs;
segMeta.probs(heuristics.segIds, 4) = heuristics.vesselScore;
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

axonssert(issorted(segGraph.edges, 2));
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

[axons, axonAgglos] = ...
    ConnectEM.buildAgglomerates( ...
        segMeta, segGraph, aggloParam.axon);

curMask = ConnectEM.checkMinLength( ...
    param, aggloParam.axon.minAggloLength, axons);
axons = axons(curMask);
axonAgglos = axonAgglos(curMask);

[dendrites, dendriteAgglos] = ...
    ConnectEM.buildAgglomerates( ...
        segMeta, segGraph, aggloParam.dend);

curMask = ConnectEM.checkMinLength( ...
    param, aggloParam.dend.minAggloLength, dendrites);
dendrites = dendrites(curMask);
dendriteAgglos = dendriteAgglos(curMask);

%% Save results
Util.log('Saving results...')
clear cur*;

curOut = struct;
curOut.info = info;
curOut.axons = axons;
curOut.dendrites = dendrites;
curOut.axonAgglos = axonAgglos;
curOut.dendriteAgglos = dendriteAgglos;

curOutFile = sprintf('%s_results.mat', datestr(now, 30));
curOutFile = fullfile(outDir, curOutFile);

Util.saveStruct(curOutFile, curOut);
keyboard
%% write skeletons
Util.log('Export skeletons:')
outputFolderSub = fullfile(outDir,'dendrites');
mkdir(outputFolderSub)
for i=1:100
    curSkel = dendrites(i);
    curSkel.write(fullfile(outputFolderSub,[sprintf('dendrite_%02d',i)]));
end

outputFolderSub = fullfile(outDir,'axons');
mkdir(outputFolderSub)
for i=1:100
    curSkel = dendrites(i);
    curSkel.write(fullfile(outputFolderSub,[sprintf('dendrite_%02d',i)]));
end
keyboard 
%% write segmentation
Util.log('Write new segmentation based on agglos')
segOut = struct;
segOut.root = fullfile(p.saveFolder, 'aggloState', ...
    [timeStamp '_segForDendAgglos'], '1');
mkdir(segOut.root)
dataType = 'uint32';
wkwInit('new',segOut.root,32, 32, dataType, 1);
segOut.backend = 'wkwrap';
agglosSorted = dendriteAgglos;
mapping = connectEM.createLookup(segMeta, agglosSorted);
Seg.Global.applyMappingToSegmentation(p, mapping, segOut);
thisBBox = [1, 1, 1; (ceil(p.bbox(:, 2) ./ 1024) .* 1024)']';
createResolutionPyramid(segOut, thisBBox, [], true);

segOut = struct;
segOut.root = fullfile(p.saveFolder, 'aggloState', ...
    [timeStamp '_segForAxonAgglos'], '1');
mkdir(segOut.root)
dataType = 'uint32';
wkwInit('new',segOut.root,32, 32, dataType, 1);
segOut.backend = 'wkwrap';
agglosSorted = axonAgglos;
mapping = connectEM.createLookup(segMeta, agglosSorted);
Seg.Global.applyMappingToSegmentation(p, mapping, segOut);
thisBBox = [1, 1, 1; (ceil(p.bbox(:, 2) ./ 1024) .* 1024)']';
createResolutionPyramid(segOut, thisBBox, [], true);
