clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
synFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v3_ax_spine_clustered.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a_ax_spine_syn_clust.mat');
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_01_spine_attachment.mat');

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

syn = load(synFile);
conn = load(connFile);

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

%%
blah = buildAxonSpineInterfaceAreas( ...
    param, graph, conn.axons, shAgglos, syn.synapses, synT);

%%