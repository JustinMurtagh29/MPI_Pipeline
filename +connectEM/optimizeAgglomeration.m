% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
%  Modified by
%   Sahil Loomba <sahil.loomba@brain.mpg.de>

%% Configuration
param = p;

gtNmlDirs = struct;
gtNmlDirs.axon = fullfile(param.saveFolder,'tracings/box-seeded/axon-ground-truth');
gtNmlDirs.dend = fullfile(param.saveFolder,'tracings/box-seeded/dendrite-ground-truth');

optimParam = struct;
optimParam.axonMinEdgeProb = 0.975;
optimParam.axonMinAxonProb = 0.5;
optimParam.axonMaxDendProb = 0.5;
optimParam.axonMaxAstroProb = 0.5;
optimParam.axonMinAggloLength = 5000;

optimParam.dendMinEdgeProb = 0.975;
optimParam.dendMinBorderSize = 0;
optimParam.dendMaxAxonProb = 0.5;
optimParam.dendMinDendProb = 0.5;
optimParam.dendMaxAstroProb = 0.5;
optimParam.dendMinAggloLength = 5000;

% Axon evaluation
axonEvalParam = struct;

% Parameters for merger detection
axonEvalParam.mergerMinDistance = 2000;
axonEvalParam.mergerRadiusInner = 1000;
axonEvalParam.mergerRadiusOuter = 10000;

axonEvalParam.nrQueriesPerSplit = 2;
axonEvalParam.tracingLenPerSplitQuery = 2000;

% NOTE(amotta): To estimate the number of merge errors let's assume that
% * the ground truth tracing is linear at all times (no branchpoints),
% * three- and four-fold chiasmata occur in equal proportion.
axonEvalParam.nrExitsPerMerger = 0.5 * (3 - 2) + 0.5 * (4 - 2);

% NOTE(amotta): Solving a chiasma / merge errors does, of course, require more
% than two queries. But these queries are "amortized" in the sense that they
% not only fix the ground truth axon, but all involved axons. So, here we
% assume two chiasma exits for the ground truth axon under consideration.
axonEvalParam.nrQueriesPerMerger = 2;

% NOTE(amotta): To first approximation, the tracing of a merge query is
% twice as long as the inner radius used for merge error detection.
axonEvalParam.tracingLenPerMergeQuery = 5*2000;


% Dendrite evaluation
dendEvalParam = struct;

% Parameters for merger detection
dendEvalParam.mergerMinDistance = 10000;
dendEvalParam.mergerRadiusInner = 5000;
dendEvalParam.mergerRadiusOuter = 20000;

dendEvalParam.nrQueriesPerSplit = 2;
dendEvalParam.tracingLenPerSplitQuery = 2000;

% NOTE(amotta): All dendritic merge errors I've seen so far were typical
% four-fold chiasmata. So, the situation is simpler than for the axons.
dendEvalParam.nrExitsPerMerger = 2;

dendEvalParam.nrQueriesPerMerger = 2;
dendEvalParam.tracingLenPerMergeQuery = 5*10000;

evalParam = struct;
evalParam.axon = axonEvalParam;
evalParam.dend = dendEvalParam;

% NOTE(amotta): The ground truth contains 6 / 69 axons and 2 / 17
% dendrites. So, let's use this info to weigh the two error classes
% accordingly.
evalParam.axonFactor = 33 / 11;
evalParam.dendFactor = 6 / 4;

parallelParam = struct;
parallelParam.batchCount = 100;
parallelParam.batchSize = 5;

info = Util.runInfo();
Util.showRunInfo(info);

%% Load data
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

%% Load ground truth skeletons
clear cur*;
gt = cell(0, 1);

curTypes = fieldnames(gtNmlDirs);
for curTypeIdx = 1:numel(curTypes)
    curType = curTypes{curTypeIdx};
    curGtNmlDir = gtNmlDirs.(curType);
    
    curNmlFiles = dir(fullfile(curGtNmlDir, '*.nml'));
    curNmlFiles(cat(1, curNmlFiles.isdir)) = [];
    curNmlFiles = {curNmlFiles.name};

    for curNmlIdx = 1:numel(curNmlFiles)
        curNmlFile = curNmlFiles{curNmlIdx};
        curNmlPath = fullfile(curGtNmlDir, curNmlFile);

        curNml = slurpNml(curNmlPath);
        curNmlNodes = NML.buildNodeTable(curNml);
        curNmlTrees = NML.buildTreeTable(curNml);

        curNmlNodes.coord = curNmlNodes.coord + 1;
        curNmlNodes.segId = Seg.Global.getSegIds(param, curNmlNodes.coord);

        curGt = struct('nodes', {}, 'segIds', {}, 'edges', {});

        for curTreeIdx = 1:height(curNmlTrees)
            curTreeId = curNmlTrees.id(curTreeIdx);
            curEdges = curNmlTrees.edges{curTreeIdx};

            curNodes = curNmlNodes.treeId == curTreeId;
            curNodes = curNmlNodes(curNodes, :);

            curEdges = cat(2, curEdges.source, curEdges.target);
           [~, curEdges] = ismember(curEdges, curNodes.id);
            assert(all(curEdges(:)));

            curEdges = reshape(curEdges, [], 2);
            curEdges = sortrows(sort(curEdges, 2));

            curGt(curTreeIdx).nodes = curNodes.coord;
            curGt(curTreeIdx).segIds = curNodes.segId;
            curGt(curTreeIdx).edges = curEdges;
            curGt(curTreeIdx).type = curType;
        end

        gt{end + 1} = curGt(:); %#ok
    end
end

gt = cat(1, gt{:});

%% Reduce multi-graph
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

%% Save shared inputs
clear cur*;

curTypes = fieldnames(gtNmlDirs);

curData = struct;
curData.param = param;
curData.segMeta = segMeta;
curData.segGraph = segGraph;
curData.gt = gt;
curData.types = curTypes(:);
curData.evalParam = evalParam;
curData.optimParam = optimParam;

curTempDir = Util.makeTempDir();
parallelParam.inputFile = fullfile(curTempDir, 'input-data.mat');

Util.log('Writing shared input data to "%s"', parallelParam.inputFile);
Util.saveStruct(parallelParam.inputFile, curData);
Util.log('Done');

%% Evaluate parameter set
clear cur*;
rng(0);

curParamNames = fieldnames(optimParam);
curNumParams = numel(curParamNames);
curLb = zeros(curNumParams, 1);
curUb = ones(curNumParams, 1);

% NOTE(amotta): Restrict edge probability thresholds to range from 50 to
% 100 %. It's highly unlikely that anything below that will make sense.
curLb(endsWith(curParamNames, 'MinEdgeProb')) = 0.5;

% NOTE(amotta): Restrict minimum agglomerate length to at most 50 µm.
curUb(endsWith(curParamNames, 'MinAggloLength')) = 50000;

% NOTE(amotta): The maximum border size exceeds unity.
curUb(endsWith(curParamNames, 'MinBorderSize')) = maxBorderSize;

curFitnessFun = @(x) fitnessParallel(parallelParam, x);
curOptions = optimoptions( ...
    'ga', 'UseVectorized', true, 'Display', 'iter', ...
    'PopulationSize', parallelParam.batchCount * parallelParam.batchSize);

[curParams, curMinCost, curExitFlags, curPop, curScores] = ga( ...
    curFitnessFun, curNumParams, [], [], [], [], ...
    curLb, curUb, [], [], curOptions);

save(fullfile(p.saveFolder,'aggloGA',[datestr(clock,30) 'optimParams.mat']),'curParams','curMinCost','curExitFlags','curPop','curScores');

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
        'memory', 48, 'cores', 4, ...
        'priority', 100, 'time', '12:00:00'};
    job = Cluster.startJob( ...
        @fitness, inputs, 'sharedInputs', sharedInputs, ...
        'numOutputs', 1, 'name', 'ga', 'cluster', clusterOpts);
    wait(job);
    
    cost = fetchOutputs(job);
    cost = cat(1, cost{:});
    cost(randIds) = cost;
end

function cost = fitness(parallelParam, inputs)
   [param, segMeta, segGraph, gt, types, evalParam, optimParam] = ...
        Util.load(parallelParam.inputFile, 'param', 'segMeta', ...
        'segGraph', 'gt', 'types', 'evalParam', 'optimParam');
    
    optimParam = fieldnames(optimParam);
    assert(isequal(numel(optimParam), size(inputs, 2)));
    optimParam = cell2struct(num2cell(inputs), optimParam, 2);
    
    % Build normalize factors
    normFactors = nan(size(types));
    for curTypeIdx = 1:numel(types)
        curType = types{curTypeIdx};
        curFactor = evalParam.(sprintf('%sFactor', curType));
        normFactors(curTypeIdx) = curFactor;
    end
    
    normFactors = normFactors / sum(normFactors);
    normFactors = cell2struct(num2cell(normFactors), types);
    
    % Evaluate parameter sets
    cost = nan(size(optimParam));
    for curParamIdx = 1:numel(optimParam)
        curOptimParam = optimParam(curParamIdx);
        curParamNames = fieldnames(curOptimParam);
        
        curEvals = struct;
        for curTypeIdx = 1:numel(types)
            curType = types{curTypeIdx};
            
            curTypeMask = reshape({gt.type}, [], 1);
            curTypeMask = strcmpi(curTypeMask, curType);
            
            % Build type-specific struct
            curTypeParam = struct;
            for curIdx = 1:numel(curParamNames)
                curOldName = curParamNames{curIdx};
                if ~startsWith(curOldName, curType); continue; end
                
                curNewName = curOldName((1 + numel(curType)):end);
                curNewName(1) = lower(curNewName(1));
                
                curTypeParam.(curNewName) = curOptimParam.(curOldName);
            end
            
            curEvalParam = evalParam.(curType);
            curEvalParam.minAggloLength = curTypeParam.minAggloLength;

            curAggloParam = rmfield(curTypeParam, 'minAggloLength');
            curAggloParam.maxVesselProb = 0.5;
            curAggloParam.maxNucleusProb = 0.5;

            curEvalT = ConnectEM.evaluateAgglomerates( ...
                param, segMeta, segGraph, gt, curEvalParam, curAggloParam);
            curEvalT.typeMask = curTypeMask;
            
            curEvals.(curType) = curEvalT;
        end
        
        curCost = 0;
        for curTypeIdx = 1:numel(types)
            curType = types{curTypeIdx};
            curEvalT = curEvals.(curType);
            curFactor = normFactors.(curType);
            
            curQueryLen = sum(curEvalT.queryLen(curEvalT.typeMask));
            curFalseLen = sum(curEvalT.pathLenCovered(~curEvalT.typeMask));
            
            curCost = curCost ...
                + curFactor * curQueryLen ...
                + (1 - curFactor) * curFalseLen;
        end
                
        cost(curParamIdx) = curCost;
    end
end

