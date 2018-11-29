% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
graphFile = '/gaba/u/hkebiri/hierarchical-agglo-graph.mat';

box = [9600, 5000, 10, 1250, 1250, 132];
box = Util.convertWebknossosToMatlabBbox(box);

param = struct;
param.bbox = box;

param.class = struct;
param.class.root = '/u/hkebiri/python/unet/unet_1fSEP_05_Noise/results/wkwPredRaw990748fSEP_05_NoiseSkeletonTracedRightDimensions';
param.class.backend = 'wkwrap';

param.seg = struct;
param.seg.root = '/tmpscratch/hkebiri/unet_1fSEP_05_Noise/segmentations_logDistTransform/h0.2_thresh_0.5_voxelSize_4435';
param.seg.backend = 'wkwrap';

param.local(1) = struct;
param.local(1).bboxSmall = box;

margin = [0; 0; 0];
% margin = param.tileBorder(:, 2);
% assert(isequal(margin, [256; 256; 128]));

info = Util.runInfo();
Util.showRunInfo(info);

%% Run per-cube agglomeration
taskArgs = arrayfun( ...
    @(local) {local.bboxSmall}, ...
    param.local, 'UniformOutput', false);
sharedArgs = {param, 'margin', margin};
sharedArgLocs = [1, 3, 4];

job = Cluster.startJob( ...
    @connectEM.Hierarchical.runBox, ...
    taskArgs, 'numOutputs', 2, ...
    'sharedInputs', sharedArgs, ...
    'sharedInputsLocation', sharedArgLocs, ...
    'cluster', {'priority', 100, 'time', '24:00:00', 'memory', 128});
Cluster.waitForJob(job);

%% Fetch outputs
out = fetchOutputs(job);
edges = cat(1, out{:, 1});
scores = cat(1, out{:, 2});
clear out;

%% Build graph
clear cur*;
[curCount, ~] = cellfun(@size, edges);
scores = repelem(scores, curCount, 1);
edges = cat(1, edges{:});

graph = table;
graph.edge = edges;
graph.score = scores;
clear edges scores;

graph = sortrows(graph, 'score', 'descend');
[~, curUni] = unique(graph.edge, 'rows', 'stable');
graph = graph(curUni, :);

Util.save(graphFile, graph);
