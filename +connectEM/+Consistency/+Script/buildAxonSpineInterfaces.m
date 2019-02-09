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

[conn, syn, connFile] = connectEM.Consistency.loadConnectome(param);

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

%% Calculate ASI areas
curBatchSize = 500;
curSharedArgs = {param, conn.axons, shAgglos};

curArgs = 1:curBatchSize:height(asiT);
curArgs = [curArgs; min(curArgs + curBatchSize - 1, height(asiT))];

curArgs = arrayfun( ...
    @(a, b) {asiT(a:b, :)}, ...
    curArgs(1, :), curArgs(2, :), ...
    'UniformOutput', false);

curJob = Cluster.startJob( ...
    @connectEM.Consistency.buildAxonSpineInterfaceAreas, curArgs, ...
    'cluster', {'priority', 100, 'memory', 24, 'cores', 2}, ...
    'numOutputs', 1, 'sharedInputs', curSharedArgs);
Cluster.waitForJob(curJob);

areas = fetchOutputs(curJob);
delete(curJob);

areas = cell2mat(areas);
asiT.area = areas;

%% Write output
clear cur*;

curOut = struct;
curOut.asiT = asiT;
curOut.info = info;

[curOutDir, curOutName] = fileparts(connFile);
curOutFile = fullfile(curOutDir, sprintf('%s_asiT.mat', curOutName));

Util.saveStruct(curOutFile, curOut);
Util.protect(curOutFile);
