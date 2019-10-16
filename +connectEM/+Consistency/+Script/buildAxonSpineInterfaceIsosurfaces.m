% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto.mat');
asiFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-linearized_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified__20190227T082543_asiT.mat');
outDir = '/tmpscratch/amotta/l4/2019-10-16-synapse-collar-isosurfaces/mat';

batchSize = 500;

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

asiT = load(asiFile, 'asiT');
asiT = asiT.asiT;

shAgglos = load(shFile, 'shAgglos');
shAgglos = shAgglos.shAgglos;

conn = connectEM.Consistency.loadConnectome(param);

%% Calculate ASI areas
clear cur;
curSharedArgs = { ...
    param, conn.axons, shAgglos, outDir};
curClusterConfig = { ...
    'priority', 100, 'memory', 24, ...
    'cores', 2, 'time', '6:00:00'};

% Restrict to interfaces with valid positions
curIds = find(not(any(isnan(asiT.pos), 2)));

curArgs = 1:batchSize:numel(curIds);
curArgs = [curArgs; min(curArgs + batchSize - 1, numel(curIds))];

curArgs = arrayfun( ...
    @(a, b) {asiT(curIds(a:b), :)}, ...
    curArgs(1, :), curArgs(2, :), ...
    'UniformOutput', false);

curJob = Cluster.startJob( ...
    @connectEM.Consistency.buildAxonSpineInterfaceIsosurfaces, ...
    curArgs, 'sharedInputs', curSharedArgs, 'cluster', curClusterConfig);
Cluster.waitForJob(curJob);
