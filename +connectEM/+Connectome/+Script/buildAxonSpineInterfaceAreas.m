clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');
% connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-linearized_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_01_spine_attachment.mat');

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;
[conn, syn] = connectEM.Connectome.load(param, connFile);

shAgglos = load(shFile, 'shAgglos');
shAgglos = shAgglos.shAgglos;

%% Loading augmented graph
graph = Graph.load(rootDir);
graph(~graph.borderIdx, :) = [];

borderAreas = fullfile(rootDir, 'globalBorder.mat');
borderAreas = load(borderAreas, 'borderArea2');
borderAreas = borderAreas.borderArea2;

graph.borderArea = borderAreas(graph.borderIdx);
clear borderAreas;

%%
synT = connectEM.Connectome.buildSynapseTable(conn, syn);

spineSynT = synT;
spineSynT.type = syn.synapses.type(spineSynT.id);
spineSynT(spineSynT.type ~= 'PrimarySpine', :) = [];

%%
asiT = connectEM.Connectome.buildAxonSpineInterfaceAreas( ...
    param, graph, conn.axons, shAgglos, syn.synapses, spineSynT);

%%
plotConfig = struct;
plotConfig.title = 'All axon-spine interface areas';
plotConfig.synIds = reshape(1:height(asiT), [], 1);

pairConfigs = connectEM.Consistency.buildPairConfigs(asiT, plotConfig);

%%
connectEM.Consistency.plotSizeHistogram( ...
    info, asiT, plotConfig, 'scale', 'log');

%%
connectEM.Consistency.plotVariabilityHistogram( ...
    info, asiT, plotConfig, pairConfigs(:));

%%
[a, b, c] = connectEM.Consistency.calculateLearnedFraction( ...
    asiT, pairConfigs(1), pairConfigs(5)) %#ok
