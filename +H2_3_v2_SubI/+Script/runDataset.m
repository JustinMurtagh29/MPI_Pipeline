% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/tmpscratch/sahilloo/data/H2_3_v2_U1_SubI/pipelineRun_mr2e_wsmrnet/';
graphFile = '/tmpscratch/sahilloo/data/H2_3_v2_U1_SubI/pipelineRun_mr2e_wsmrnet/29Nov2018_agglomeration/graph.mat';

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

% add correspondence scores
corrFile = fullfile(param.saveFolder, 'mapping.mat');
corr = load(corrFile, 'correspondences');
corr = corr.correspondences;

borderIdx = [
        reshape(1:size(edges, 1), [], 1); % for border-based edges
        nan(size(corr, 1), 1)];           % for correspondences
    
edges = [edges; corr];
scores = [scores; ones(size(corr, 1), 1)];

% sort edges and probabilities
[edges, rows] = sortrows(edges);
scores = scores(rows);
borderIdx = borderIdx(rows);

% build output
graph = table;
graph.edges = edges;
graph.scores = scores;
graph.borderIdx = borderIdx;

% save output
graph = sortrows(graph, 'scores', 'descend');
[~, curUni] = unique(graph.edges, 'rows', 'stable');
graph = graph(curUni, :);

Util.save(graphFile, graph, info)




