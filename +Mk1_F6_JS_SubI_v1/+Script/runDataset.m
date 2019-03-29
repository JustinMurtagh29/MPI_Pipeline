% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
datasetName = 'Mk1_F6_JS_SubI_v1';
rootDir = fullfile('/tmpscratch/sahilloo/data',datasetName,'pipelineRun_mr2e_wsmrnet');
graphFile = fullfile(rootDir,'agglomeration','graph.mat');

param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

margin = param.tileBorder(:, 2);
assert(isequal(margin, [256; 256; 128]));

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
    'cluster', {'priority', 100, 'time', '24:00:00', 'memory', 48});
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
