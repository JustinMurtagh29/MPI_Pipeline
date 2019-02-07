% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto.mat');

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, syn, connFile] = ...
    connectEM.Consistency.loadConnectome(param);

% Loading spine head agglomerates
shAgglos = load(shFile, 'shAgglos');
shAgglos = shAgglos.shAgglos;

borderMeta = fullfile(rootDir, 'globalBorder.mat');
borderMeta = load(borderMeta, 'borderArea2', 'borderSize', 'borderCoM');
borderMeta = structfun(@double, borderMeta, 'UniformOutput', false);

% Loading augmented graph
graph = Graph.load(rootDir);
graph = graph(graph.borderIdx ~= 0, :);
graph.borderArea = borderMeta.borderArea2(graph.borderIdx);

%% Build axon-spine interface areas
asiT = ...
    connectEM.Connectome.buildAxonSpineInterfaces( ...
        param, graph, shAgglos, conn, syn, 'addBorderIdsVar', true);

%% Calculate ASI positions
weightedMean = @(w, v) ...
    sum((w / sum(w, 1)) .* v, 1);

asiT.pos = cellfun( ...
    @(ids) weightedMean( ...
        borderMeta.borderSize(ids), ...
        borderMeta.borderCoM(ids, :)), ...
	asiT.borderIds, 'UniformOutput', false);
asiT.pos = round(cell2mat(asiT.pos));

%% Calculate areas
curBatch = 500;
curSharedArgs = {param, conn.axons, shAgglos};

curArgs = 1:curBatch:height(asiT);
curArgs = [curArgs; min(curArgs + curBatch - 1, height(asiT))];

curArgs = arrayfun( ...
    @(a, b) {asiT(a:b, :)}, ...
    curArgs(1, :), curArgs(2, :), ...
    'UniformOutput', false);

job = Cluster.startJob( ...
    @connectEM.Consistency.axonSpineInterfaceArea, curArgs, ...
    'cluster', {'priority', 100, 'memory', 24, 'cores', 2}, ...
    'numOutputs', 1, 'sharedInputs', curSharedArgs);
Cluster.waitForJob(job);

areas = fetchOutputs(job);
areas = cell2mat(areas);
