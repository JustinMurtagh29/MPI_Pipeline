rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-linearized_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

graphFile = fullfile(rootDir, 'graph.mat');
graph = load(graphFile, 'edges', 'borderIdx');

% NOTE(amotta): Use border idx zero instead of nan
graph.borderIdx(isnan(graph.borderIdx)) = 0;

save('/gaba/u/yyener/astrocyte/synapses/graph.mat', 'graph')
