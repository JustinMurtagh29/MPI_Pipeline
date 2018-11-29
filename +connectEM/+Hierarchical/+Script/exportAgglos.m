% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/amotta/data/2018-11-29-msem-hierarchical-agglomeration-graph';
outFile = '/home/amotta/Desktop/test.nml';

% For export to webKnossos
datasetName = '2017-05-31_mSEM_SegU-Net-Check-Mergers_3';
voxelSize = [4, 4, 35];

% Agglomerate down to this score treshold
minScore = 0;

info = Util.runInfo();

%% Loading data
meta = load(fullfile(rootDir, 'segmentMeta.mat'));
maxSegId = meta.maxSegId;
points = transpose(meta.point);
clear meta;

graph = load(fullfile(rootDir, 'graph.mat'));
graph = graph.graph;

%% Build agglomerates and export skeleton
mergeEdges = graph.edge(graph.score > minScore, :);
[~, agglos] = Graph.buildConnectedComponents(maxSegId, mergeEdges);

% HACK(amotta): Convert graph table into historical graph table;
graphS = struct;
graphS.edges = graph.edge;

skel = Skeleton.fromAgglo(graphS, points, agglos);
skel = skel.setParams(datasetName, voxelSize, [0, 0, 0]);

description = sprintf('%s (%s)', info.filename, info.git_repos{1}.hash);
skel = skel.setDescription(description);

skel.write(outFile);
