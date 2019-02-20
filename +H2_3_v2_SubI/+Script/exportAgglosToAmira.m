% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
% Modified by
%   Sahil Loomba <sahil.loomba@brain.mpg.de>
clear;

%% Configuration
rootDir = '/tmpscratch/sahilloo/data/H2_3_v2_U1_SubI/pipelineRun_mr2e_wsmrnet/';

% For export to webKnossos
load(fullfile(rootDir,'/allParameter.mat'))
datasetName = 'H2_3_v2_U1_SubI_mr2e_wsmrnet';
voxelSize = [11.24, 11.24, 28];

% Agglomerate down to this score treshold
minScore = 0.1;

info = Util.runInfo();

%% Loading data
meta = load(fullfile(rootDir, 'segmentMeta.mat'));
maxSegId = meta.maxSegId;
points = transpose(meta.point);
clear meta;
Util.log('loading graph...')
graph = load(fullfile(rootDir, '29Nov2018_agglomeration/graph.mat'));
graph = graph.graph;

Util.log('Build agglomerates and export skeleton')
mergeEdges = graph.edge(graph.score > minScore, :);
[~, agglos] = Graph.buildConnectedComponents(maxSegId, mergeEdges);

len = cellfun(@(x) length(x),agglos);
[~,idxSort] = sort(len,'descend');
agglos = agglos(idxSort);

% HACK(amotta): Convert graph table into historical graph table;
graphS = struct;
graphS.edges = graph.edge;

Util.log('Now writing out nmls...')
outDir = fullfile(rootDir,'29Nov2018_agglomeration', ['score_', num2str(minScore) '_isos']);
mkdir(outDir)

agglosOut = agglos(2:300);

Visualization.exportAgglosToAmira(p, agglosOut, outDir, 'reduce',0.05,'smoothSizeHalf',4,'smoothWidth',8)
