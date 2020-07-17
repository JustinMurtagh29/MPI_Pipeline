% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
% Modified by
%   Sahil Loomba <sahil.loomba@brain.mpg.de>

% This module was developed by AM and modifed here for NHP dataset
% Note: replace astrocytes class with glia
%       heuristics: update nans error

%{
clear;

%% Configuration
rootDir = '/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet/';
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

nmlDir = fullfile(param.saveFolder, 'tracings', 'spine-attachment');
dendFile = fullfile(param.saveFolder, 'aggloMat', '20191227T134319_ga_20191224T235355optimParams_agglomeration/20191227T220548_results.mat');

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

%% Loading node and edge data
Util.log('Loading graph');
graph = fullfile(param.saveFolder, 'graph.mat');
graph = load(graph, 'edges', 'prob');

Util.log('Finding graph neighbors');
graph = Graph.addNeighbours(graph);

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
%}
for curIdx = 1:numel(nmlFiles)
    curNmlFile = nmlFiles{curIdx};
    spineT = mSEM.Spine.Attachment.loadNml(param, curNmlFile);
end

% Generate consecutive spine IDs
[~, ~, spineT.id] = unique(spineT.id);

%% Build spine head agglomerates
clear cur*;

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
gtT(cellfun(@isempty, gtT.trunkAgglo), :) = [];

% NOTE(amotta): Let's also ignore the ground truth tracings where none of
% the trunk nodes was part of a dendrite agglomerate.
curMask = cellfun(@isempty, gtT.trunkIds);
warning('Ignoring %d spines that did not reach a trunk', sum(curMask));
gtT(curMask, :) = [];

%% Evaluate parameter set
clear cur*;
rng(0);

curData = struct;
curData.graph = graph;
curData.dendLUT = dendLUT(:);
curData.segMeta = segMeta;

curParamNames = fieldnames(optimParam);
curParamVals = cell2mat(struct2cell(optimParam));
curParamVals = reshape(curParamVals, 1, []);

[curCost, curUnattached, curWronglyAttached] = ...
    fitness(gtT, curData, curParamNames, curParamVals);

%% Core
function [cost, numUnattached, numWronglyAttached] = ...
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
    
    attachData.exclLUT = ...
        (curProbs(:, curClasses == 'axon') > attachParams.maxAxonProb) ...
      | (curProbs(:, curClasses == 'dendrite') < attachParams.minDendProb) ...
      | (curProbs(:, curClasses == 'glia') > attachParams.maxAstroProb) ...
      | (curProbs(:, curClasses == 'vessel') > 0.5) ...
      | (curProbs(:, curClasses == 'nucleus') > 0.5);
  
   [~, attachedTo] = ...
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
