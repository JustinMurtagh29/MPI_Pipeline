% This script is a modified version of
% +connectEM/+Consistency/+Script/buildAxonSpineInterfaces.m
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet/';
connFile = fullfile(rootDir, 'connectome', 'Connectome_20191227T220548-results_20191227T220548-results-auto-spines-v3_SynapseAgglomerates--20191227T220548-results--20191227T220548-results-auto-spines-v3--v1.mat');

runId = datestr(now, 30);

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

conn = load(connFile);
% NOTE(amotta): We haven't defined any axon classes yet.
conn.axonMeta.axonClass(:) = categorical({'Unknown'});

syn = load(conn.info.param.synFile);
% NOTE(amotta): We do not yet distinguish between so-called primary and
% secondary spines. So, let's call our categories just "Spine" and "Shaft".
syn.synapses.type = categorical( ...
    1 + syn.isSpineSyn, [1, 2], {'Shaft', 'Spine'});

shFile = syn.info.param.spineHeadFile;

% Loading spine head agglomerates
shAgglos = load(shFile, 'shAgglos');
shAgglos = shAgglos.shAgglos;

borderMeta = fullfile(rootDir, 'globalBorder.mat');
borderMeta = load(borderMeta, 'borderSize', 'borderCoM');
borderMeta = structfun(@double, borderMeta, 'UniformOutput', false);

% Loading augmented graph
graph = Graph.load(rootDir);
graph = graph(graph.borderIdx ~= 0, :);

%% Build axon-spine interfaces (ASIs)
asiT = ...
    connectEM.Connectome.buildAxonSpineInterfaces( ...
        param, graph, shAgglos, conn, syn, 'addBorderIdsVar', true);

asiT.axonClass = conn.axonMeta.axonClass(asiT.preAggloId);
asiT.targetClass = conn.denMeta.targetClass(asiT.postAggloId);

%% Calculate ASI positions
weightedMean = @(w, v) ...
    sum((w / sum(w, 1)) .* v, 1);

asiT.pos = cellfun( ...
    @(ids) weightedMean( ...
        borderMeta.borderSize(ids), ...
        borderMeta.borderCoM(ids, :)), ...
	asiT.borderIds, 'UniformOutput', false);
asiT.pos = round(cell2mat(asiT.pos));

% NOTE(amotta): The synapse agglomeration happened independently of the
% neurite reconstruction. In very rare cases (55 times at the time of
% writing) it can thus happen that there does not exist any border
% between a axon and spine head agglomerates.
%   If the set of border IDs is empty, the empty sum in the weighted
% mean above returns all zeros. Let's set the position of these alien
% interfaces to not-a-number instead.
asiT.pos(cellfun(@isempty, asiT.borderIds), :) = nan;

%% Calculate ASI areas
curBatchSize = 500;
curSharedArgs = {param, conn.axons, shAgglos};

% Restrict to interfaces with valid positions
curIds = find(not(any(isnan(asiT.pos), 2)));

curArgs = 1:curBatchSize:numel(curIds);
curArgs = [curArgs; min(curArgs + curBatchSize - 1, numel(curIds))];

curArgs = arrayfun( ...
    @(a, b) {asiT(curIds(a:b), :)}, ...
    curArgs(1, :), curArgs(2, :), ...
    'UniformOutput', false);

curJob = Cluster.startJob( ...
    @connectEM.Consistency.buildAxonSpineInterfaceAreas, curArgs, ...
    'cluster', {'priority', 100, 'memory', 24, 'cores', 2, 'time', '6:00:00'}, ...
    'numOutputs', 1, 'sharedInputs', curSharedArgs);
Cluster.waitForJob(curJob);

areas = fetchOutputs(curJob);
areas = cell2mat(areas);
delete(curJob);

asiT.area = nan(height(asiT), 1);
asiT.area(curIds) = areas;

%% Write output
clear cur*;

curOut = struct;
curOut.asiT = asiT;
curOut.info = info;

[curOutDir, curOutName] = fileparts(connFile);
curOutFile = sprintf('%s__%s_asiT.mat', curOutName, runId);
curOutFile = fullfile(curOutDir, curOutFile);

Util.saveStruct(curOutFile, curOut);
Util.protect(curOutFile);
