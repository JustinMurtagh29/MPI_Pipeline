% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/tmpscratch/sahilloo/data/H2_3_v2_U1_SubI/pipelineRun_mr2e_wsmrnet/';
outFile = '/tmpscratch/sahilloo/data/H2_3_v2_U1_SubI/pipelineRun_mr2e_wsmrnet/29Nov2018_agglomeration/test2.nml';

% For export to webKnossos
load(fullfile(rootDir,'/allParameter.mat'))
datasetName = 'H2_3_v2_U1_SubI_mr2e_wsmrnet';
voxelSize = [11.24, 11.24, 28];

% Agglomerate down to this score treshold
minScore = -0.90;

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

% look at random 100 agglos
idxOut = randperm(size(agglos,1),100);
agglosOut = agglos(idxOut);
% HACK(amotta): Convert graph table into historical graph table;
graphS = struct;
graphS.edges = graph.edge;

skel = Skeleton.fromAgglo(graphS, points, agglosOut);
skel = skel.setParams(datasetName, voxelSize, [0, 0, 0]);

description = sprintf('%s (%s)', info.filename, info.git_repos{1}.hash);
skel = skel.setDescription(description);

skel.write(outFile);

