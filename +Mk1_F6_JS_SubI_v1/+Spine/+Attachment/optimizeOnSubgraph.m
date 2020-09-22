% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
% Modified by
%   Sahil Loomba <sahil.loomba@brain.mpg.de>

% This module was developed by AM and modifed here for NHP dataset
% Note: replace astrocytes class with glia
%       heuristics: update nans error

clear;

%% Configuration
rootDir = '/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet/';
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

nmlDir = fullfile(param.saveFolder, 'tracings', 'spine-attachment');
dendFile = fullfile(param.saveFolder, 'aggloMat', '20191227T134319_ga_20191224T235355optimParams_agglomeration/20191227T220548_results.mat');
outFile =  fullfile(param.saveFolder,'spineAttachment',['optimizeOnSubgraph-connectem-spines-' datestr(clock, 30) '.mat']);

boxPadNm = 10e3; %5000
maxNumSteps = 20;
minEdgeProb = 0.5;

optimParam = struct;
% NOTE(amotta): These are the mean values over the final population of 100 in
% which all parameter sets reached a minimal cost of 9. Let's check whether the
% mean performs as the individual parameter sets.
optimParam.minEdgeProb = 0.0041;
optimParam.maxNumSteps = 12.20;
optimParam.maxAxonProb = 0.9510;
optimParam.minDendProb = 0.1659;
optimParam.maxAstroProb = 0.5434;

info = Util.runInfo();
Util.showRunInfo(info);

% NOTE(amotta): This field is assumed to exist by some of the functions
% that were originally written for the analysis of mSEM data (e.g.,
% mSEM.Spine.Attachment.loadNml)
param.mag = [1, 1, 1];

maxSegId = Seg.Global.getMaxSegId(param);

%% Loading dendrites
dendSuperAgglos = ...
    Util.load(dendFile, 'dendrites');
dendAgglos = arrayfun( ...
    @(s) s.segIds, dendSuperAgglos, ...
    'UniformOutput', false);

dendLUT = Agglo.buildLUT(maxSegId, dendAgglos);

% COPYNPASTE(amotta): The following code was copied from
% +L23/+ConnectEM/optimizeAgglomeration.m
segMeta = struct;
segMeta.id = reshape(1:maxSegId, [], 1);

Util.log('Loading segment types');
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
segMeta.probs(:, 4) = zeros(maxSegId,1); % all non vessel seg have score 0, not NaN
segMeta.probs(heuristics.segIds, 4) = heuristics.vesselScore;
segMeta.probs(:, 5) = zeros(maxSegId,1); % all non nuclei seg have score 0, not NaN
segMeta.probs(heuristics.segIds, 5) = heuristics.nucleiScore;

segMeta.classes = categorical([ ...
    cellstr(types.classes), {'vessel', 'nucleus'}]);

clear types

%% Loading ground truth tracings
clear cur*;

nmlFiles = dir(fullfile(nmlDir, '*.nml'));
nmlFiles(cat(1, nmlFiles.isdir)) = [];
nmlFiles = fullfile(nmlDir, {nmlFiles.name});

spineT = table();
cubeIds = [];
for curIdx = 1:numel(nmlFiles)
    curNmlFile = nmlFiles{curIdx};
    [tempSpineT, boxMag1] = mSEM.Spine.Attachment.loadNml(param, curNmlFile);
    if ~isempty(spineT)
        tempSpineT.id = max(spineT.id) + tempSpineT.id; % append to spineT
    end
    spineT = cat(1, spineT, tempSpineT);

    % padded box
    box = ceil(boxMag1 ./ param.mag(:));
    boxPadVx = ceil(boxPadNm ./ param.raw.voxelSize);
    boxPadded = box + [-1, +1] .* boxPadVx(:);

    % cubes to load graph from
    cubeIds = cat(1, cubeIds, unique(Util.boxToCubeIds(param, boxPadded)));

    clear tempSpineT, boxMag1, boxPadded
end

% Generate consecutive spine IDs
[~, ~, spineT.id] = unique(spineT.id);

%% Build spine head agglomerates
gtT = table;
gtT.id = reshape(1:max(spineT.id), [], 1);

curTemp = spineT( ...
    spineT.type == 'Head', :);
gtT.shAgglo = accumarray( ...
    curTemp.id, curTemp.segId, ...
    [], @(ids) {unique(ids(:))}, ...
    {zeros(0, 1, 'like', curTemp.segId)});

curTemp = spineT( ...
    spineT.type == 'Trunk', :);
gtT.trunkAgglo = accumarray( ...
    curTemp.id, curTemp.segId, [], ...
    @(ids) {reshape(setdiff(ids, 0), [], 1)}, ...
    {zeros(0, 1, 'like', curTemp.segId)});
gtT.trunkIds = cellfun( ...
    @(ids) reshape(setdiff(dendLUT(ids), 0), [], 1), ...
    gtT.trunkAgglo, 'UniformOutput', false);

% NOTE(amotta): For now, let's ignore ground truth tracings in which even
% manual attachment was not possible (for me, at least).
curMask = cellfun(@isempty, gtT.trunkAgglo);
warning('Ignoring %d spines where no trunk was found manually', sum(curMask));
gtT(curMask,:) = [];

% NOTE(amotta): Let's also ignore the ground truth tracings where none of
% the trunk nodes was part of a dendrite agglomerate.
curMask = cellfun(@isempty, gtT.trunkIds);
warning('Ignoring %d spines where trunk segId based dendrite agglo was not found', sum(curMask));
gtT(curMask, :) = [];

%% Load graph around box
Util.log('Loading subgraph')
edges = cell(numel(cubeIds), 1);
probs = cell(numel(cubeIds), 1);

for curIdx = 1:numel(cubeIds)
    curCubeId = cubeIds(curIdx);
    curCubeRoot = param.local(curCubeId).saveFolder;

    curProbFile = fullfile(curCubeRoot, 'neuriteContinuityProb.mat');
    curEdgeFile = fullfile(curCubeRoot, 'edges.mat');

    curProbs = Util.load(curProbFile, 'prob');
    curEdges = Util.load(curEdgeFile, 'edges');

    edges{curIdx} = curEdges;
    probs{curIdx} = curProbs;
end

graph = table;
edges = vertcat(edges{:});
graph.edges = edges;
probs = vertcat(probs{:});
graph.prob = probs;
clear edges probs;

Util.log('Post-processing graph');
uniSegIds = unique(graph.edges(:));
[~, graph.edges] = ismember(graph.edges, uniSegIds);

graph = Graph.reduce(graph);
%graph = graph(graph.prob > minEdgeProb, :);

Util.log('Adding neighbors');
graph = table2struct(graph, 'ToScalar', true);
graph = Graph.addNeighbours(graph);

% Remove shAgglos not in the graph neighbourhood
shAgglos = cellfun(@(x) x(1), gtT.shAgglo);
[curMask, shAgglos] = ismember(shAgglos, uint64(uniSegIds));
shAgglos = num2cell(shAgglos);
gtT.shAgglo = shAgglos;
warning(sprintf('Removed %d shAgglos not found in sub-graph',sum(~curMask)))

%% Evaluate parameter set
Util.log('Evaluation parameter set:')
rng(0);

curData = struct;
curData.graph = graph;
curData.dendLUT = false(numel(uniSegIds), 1); %dendLUT(:);
curData.segMeta = segMeta;
curData.exclLUT = false(numel(uniSegIds), 1);

curParamNames = fieldnames(optimParam);
curParamVals = cell2mat(struct2cell(optimParam));
curParamVals = reshape(curParamVals, 1, []);

[curCost, curUnattached, curWronglyAttached, edges, attached] = ...
    fitness(gtT, curData, curParamNames, curParamVals);

% Convert back to global segment IDs
for i = 1:numel(edges)
    curEdges = edges{i};
    curMask = curEdges ~= 0;
    curEdges(curMask) = uniSegIds(curEdges(curMask));
    edges{i} = curEdges;
end

%% Build NML file
Util.log('Writing nmls:')
rng(0);

curNumSpineHeads = 100;
curNumSpineHeads = min(curNumSpineHeads, numel(edges));
curNumDigits = ceil(log10(1 + curNumSpineHeads));

curRandIds = randperm(numel(edges));
curRandIds = curRandIds(1:curNumSpineHeads);

if ~isempty(outFile)
    segmentMeta = load([param.saveFolder 'segmentMeta.mat'], 'voxelCount', 'point', 'maxSegId');
    segmentMeta.point = segmentMeta.point';

    skel = skeleton();
    skel = Skeleton.setParams4Pipeline(skel, param);
    skel = Skeleton.setDescriptionFromRunInfo(skel, info);

    for curIdx = 1:curNumSpineHeads
        curId = curRandIds(curIdx);
        
        curSegIds = setdiff(transpose(edges{curId}), 0, 'stable');
        curSegIds = union(uniSegIds(shAgglos{curId}), curSegIds, 'stable');
        
        %curSegPos = mSEM.Seg.getPointsBySegIds(param, curSegIds);
        curSegPos = segmentMeta.point(curSegIds,:);        
        %curSegPos = transpose(curSegPos);

       [~, curEdges] = ismember(edges{curId}, curSegIds);
        curEdges = curEdges(all(curEdges, 2), :);

        curComments = repelem({''}, numel(curSegIds));
        curComments{1} = 'Spine head';

        curName = sprintf( ...
            '%0*d. Spine head %d', ...
            curNumDigits, curIdx, curId);
        skel = skel.addTree( ...
            curName, curSegPos, curEdges, ...
            [], [], curComments);
    end
    
    skel.write(outFile);
 
end

%% Core
function [cost, numUnattached, numWronglyAttached, edges, attachedTo] = ...
        fitness(gtT, data, paramNames, paramVals)
    assert(isequal(numel(paramNames), numel(paramVals)));
    attachParams = cell2struct(num2cell(paramVals), paramNames, 2);

    % HACKHACKHACK(amotta): Make sure it's an integer
    attachParams.maxNumSteps = round(attachParams.maxNumSteps);
    
    attachData = struct;
    attachData.graph = data.graph;
    attachData.dendLUT = data.dendLUT;
    
    curProbs = data.segMeta.probs;
    curClasses = data.segMeta.classes;
    
    attachData.exclLUT = data.exclLUT;
   
    [edges, attachedTo] = ...
        L4.Spine.Head.attachCore( ...
            attachParams, attachData, gtT.shAgglo);

    isAttached = attachedTo > 0;
    isCorrect = cellfun(@ismember, num2cell(attachedTo), gtT.trunkIds);
    assert(not(any(isCorrect(not(isAttached)))));

    % NOTE(amotta): We punish wrongly attached spine heads twice as much as
    % non-attached spine heads. The factor is, of course, arbitrary.
    numUnattached = sum(not(isAttached));
    numWronglyAttached = sum(not(isCorrect(isAttached)));
    cost = numUnattached + 2 * numWronglyAttached;
end
