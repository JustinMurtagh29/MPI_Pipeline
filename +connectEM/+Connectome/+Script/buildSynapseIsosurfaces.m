% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');
outDir = '/tmpscratch/amotta/l4/2020-02-06-synapse-collar-isosurfaces/mat';

batchSize = 500;
minSynPre = 10;

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, syn, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);
[conn, axonClasses] = ...
    connectEM.Connectome.prepareForSpecificityAnalysis( ...
        conn, axonClasses, 'minSynPre', minSynPre);
    
%% Prepare synapse table
curSynPos = connectEM.Synapse.calculatePositions(param, syn, 'border');
curSynPos = round(curSynPos);

synT = connectEM.Connectome.buildSynapseTable(conn, syn);
synT.pos = curSynPos(synT.id, :);

%% Calculate ASI areas
clear cur;
curSharedArgs = { ...
    param, conn.axons, ...
    conn.dendrites, outDir};
curClusterConfig = { ...
    'priority', 100, 'memory', 24, ...
    'cores', 2, 'time', '6:00:00'};

curArgs = arrayfun( ...
    @(i) {synT(i, :)}, 1:height(synT), ...
    'UniformOutput', false);

curJob = Cluster.startJob( ...
    @connectEM.Connectome.buildSynapseIsosurfaces, ...
    curArgs, 'sharedInputs', curSharedArgs, ...
    'cluster', curClusterConfig, ...
    'taskGroupSize', batchSize);
Cluster.waitForJob(curJob);
