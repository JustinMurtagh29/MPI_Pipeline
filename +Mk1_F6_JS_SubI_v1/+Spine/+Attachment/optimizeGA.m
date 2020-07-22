% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
% Modified by
%   Sahil Loomba <sahil.loomba@brain.mpg.de>

% This module was developed by AM and modifed here for NHP dataset
% Note: replace astrocytes class with glia
%       heuristics: update nans error
% use GA to optimize the attachment parameters
%{
clear;

%% Configuration
rootDir = '/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet/';
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

nmlDir = fullfile(param.saveFolder, 'tracings', 'spine-attachment');
dendFile = fullfile(param.saveFolder, 'aggloMat', '20191227T134319_ga_20191224T235355optimParams_agglomeration/20191227T220548_results.mat');

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

parallelParam = struct;
parallelParam.batchCount = 100;
parallelParam.batchSize = 5;

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

graph = struct2table(graph);
graph = Graph.reduce(graph);
graph = graph(graph.prob > minEdgeProb, :);

Util.log('Finding graph neighbors');
graph = table2struct(graph, 'ToScalar', true);
%graph = Graph.addNeighbours(graph);

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
nmlFiles = dir(fullfile(nmlDir, '*.nml'));
nmlFiles(cat(1, nmlFiles.isdir)) = [];
nmlFiles = fullfile(nmlDir, {nmlFiles.name});

spineT = table();
for curIdx = 1:numel(nmlFiles)
    curNmlFile = nmlFiles{curIdx};
    tempSpineT = mSEM.Spine.Attachment.loadNml(param, curNmlFile);
    if ~isempty(spineT)
        tempSpineT.id = max(spineT.id) + tempSpineT.id; % append to spineT
    end
    spineT = cat(1, spineT, tempSpineT);
    clear tempSpineT
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

%% Save shared inputs
curData = struct;
curData.gtT = gtT;
curData.graph = graph;
curData.dendLUT = dendLUT(:);
curData.segMeta = segMeta;
curData.optimParam = optimParam;

%parallelParam.graph = graph; % graph too big to save and load from disk
curTempDir = Util.makeTempDir();
parallelParam.inputFile = fullfile(curTempDir, 'input-data.mat');

Util.log('Writing shared input data to "%s"', parallelParam.inputFile);
Util.saveStruct(parallelParam.inputFile, curData);
Util.log('Done');
%}
%% Evaluate parameter set
rng(0);

curParamNames = fieldnames(optimParam);
curNumParams = numel(curParamNames);
curLb = zeros(curNumParams, 1);
curUb = ones(curNumParams, 1);

% It's highly unlikely that anything below that will make sense.
%curLb(endsWith(curParamNames, 'minEdgeProb')) = 0.5;

% allow max steps to 15
curUb(endsWith(curParamNames, 'maxNumSteps')) = maxNumSteps;

curFitnessFun = @(x) fitnessParallel(parallelParam, x);
curOptions = optimoptions( ...
    'ga', 'UseVectorized', true, 'Display', 'iter', ...
    'PopulationSize', parallelParam.batchCount * parallelParam.batchSize);

[curParams, curMinCost, curExitFlags, curPop, curScores] = ga( ...
    curFitnessFun, curNumParams, [], [], [], [], ...
    curLb, curUb, [], [], curOptions);

outfile = fullfile(param.saveFolder,'spineAttachment',[datestr(clock,30) 'optimParam.mat']);
save(outfile, 'curParams','curMinCost','curExitFlags','curPop','curScores');
Util.log(sprintf('Finished saving optimized curParams to file \n %s', outfile))

%% Utility
function cost = fitnessParallel(parallelParam, inputs)
    % NOTE(amotta): Randomize order of parameter sets, so that run times
    % average out over tasks.
    randIds = randperm(size(inputs, 1));
    inputs = inputs(randIds, :);

    inputs = mat2cell(inputs, repelem(parallelParam.batchSize, ...
        parallelParam.batchCount), size(inputs, 2));
    inputs = arrayfun(@(x) x, inputs, 'UniformOutput', false);
    
    sharedInputs = {parallelParam};
    
    clusterOpts = { ...
        'memory', 72, 'cores', 6, ...
        'priority', 100, 'time', '12:00:00'};
    job = Cluster.startJob( ...
        @fitness, inputs, 'sharedInputs', sharedInputs, ...
        'numOutputs', 1, 'name', 'ga', 'cluster', clusterOpts);
    wait(job);
    
    cost = fetchOutputs(job);
    cost = cat(1, cost{:});
    cost(randIds) = cost;
end

%% Core
function cost = fitness(parallelParam, inputs)
    % graph is too big to store and read
    %graph = parallelParam.graph;

    [gtT, dendLUT, graph, segMeta, optimParam] = Util.load(parallelParam.inputFile, ...
            'gtT', 'dendLUT', 'graph', 'segMeta', 'optimParam');
    graph = Graph.addNeighbours(graph);
 
    optimParam = fieldnames(optimParam);
    assert(isequal(numel(optimParam), size(inputs, 2)));
    optimParam = cell2struct(num2cell(inputs), optimParam, 2);

    % evaluate parameter sets
    cost = nan(size(optimParam));
    for curParamIdx = 1:numel(optimParam)
        curOptimParam = optimParam(curParamIdx);

        % HACKHACKHACK(amotta): Make sure it's an integer
        curOptimParam.maxNumSteps = round(curOptimParam.maxNumSteps);
        
        attachData = struct;
        attachData.graph = graph;
        attachData.dendLUT = dendLUT;
        
        curProbs = segMeta.probs;
        curClasses = segMeta.classes;
        
        attachData.exclLUT = ...
            (curProbs(:, curClasses == 'axon') > curOptimParam.maxAxonProb) ...
          | (curProbs(:, curClasses == 'dendrite') < curOptimParam.minDendProb) ...
          | (curProbs(:, curClasses == 'glia') > curOptimParam.maxAstroProb) ...
          | (curProbs(:, curClasses == 'vessel') > 0.5) ...
          | (curProbs(:, curClasses == 'nucleus') > 0.5);
      
        [~, attachedTo] = ...
            L4.Spine.Head.attachCore( ...
                curOptimParam, attachData, gtT.shAgglo);
    
        isAttached = attachedTo > 0;
        isCorrect = cellfun(@ismember, num2cell(attachedTo), gtT.trunkIds);
        assert(not(any(isCorrect(not(isAttached)))));
    
        % NOTE(amotta): We punish wrongly attached spine heads twice as much as
        % non-attached spine heads. The factor is, of course, arbitrary.
        numUnattached = sum(not(isAttached));
        numWronglyAttached = sum(not(isCorrect(isAttached)));
        cost(curParamIdx) = numUnattached + 2 * numWronglyAttached;
    end
end
