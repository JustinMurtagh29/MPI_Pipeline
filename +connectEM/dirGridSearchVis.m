%{
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
%}

% Now everything including old and some new metrics, see:
% https://mhlablog.net/2017/05/16/focusem-path-length-leverage/

functionH = @connectEM.moreMetrics;
inputArguments = arrayfun(@(x,y){x,y{1}}, runs, results', 'uni', 0);
cluster = Cluster.getCluster( ...
    '-pe openmp 1', ...
    '-p -500', ...
    '-l h_vmem=40G', ... 
    '-l s_rt=99:50:00', ...
    '-l h_rt=100:00:00');
job = Cluster.startJob( functionH, inputArguments, ...
    'name', 'moreMetrics', ...
    'sharedInputs', {workingFolder segmentMeta}, ...
    'sharedInputsLocation', [1 4], ...
    'numOutputs', 1, ...
    'cluster', cluster);

Cluster.waitForJob(job);
metrics = fetchOutputs(job);

