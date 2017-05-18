% Load needed information
graph = load('/gaba/u/mberning/results/pipeline/20170217_ROI/graphNew.mat', 'edges', 'prob', 'borderIdx');
[graph.neighbours, neighboursIdx] = Graph.edges2Neighbors(graph.edges);
graph.neighProb = cellfun(@(x)graph.prob(x), neighboursIdx, 'uni', 0);
graph.neighBorderIdx = cellfun(@(x)graph.borderIdx(x), neighboursIdx, 'uni', 0);
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
segmentMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta.mat');
segmentMeta.point = segmentMeta.point';
segmentMeta = connectEM.addSegmentClassInformation(p, segmentMeta);

% Where to look, what is there
workingFolder = '/gaba/scratch/mberning/aggloGridSearch2/';
load([workingFolder 'parameters.mat']);
runs = dir([workingFolder 'search01_*']);
for i=1:length(runs)
    results{i} = dir([workingFolder runs(i).name filesep '*.mat']);
end

% Now everything including old and some new metrics, see:
% https://mhlablog.net/2017/05/16/focusem-path-length-leverage/

% Metrics on all runs
tic;
for i=1:length(results)
    metrics(i).nrRuns = numel(results{i});
    load([workingFolder runs(i).name filesep results{i}(metrics.nrRuns).name], 'axonsNew');
    metrics(i).y = connectEM.evaluateAggloMetaMeta(graph, axonsNew, [], [runs(i).name '_' metrics.nrRuns], segmentMeta);
    metrics(i).pathLength  = connectEM.getPathLengthFromAgglomeration(axonsNew, segmentMeta.point);
    metrics(i).nrLargeAgglos = sum(metrics(i).pathLength > 5);
    metrics(i).pathLengthLargeAgglos = sum(metrics(i).pathLength(metrics(i).pathLength > 5));
    metrics(i).yLarge = connectEM.evaluateAggloMetaMeta(graph, axonsNew(metrics(i).pathLength > 5), [], [runs(i).name '_' metrics.nrRuns '_large'], segmentMeta);
    Util.progressBar(i, length(results));
end

