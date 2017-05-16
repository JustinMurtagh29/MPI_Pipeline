% Load Graph
graph = load('/gaba/u/mberning/results/pipeline/20170217_ROI/graphNew.mat', 'edges', 'prob', 'borderIdx');
[graph.neighbours, neighboursIdx] = Graph.edges2Neighbors(graph.edges);
graph.neighProb = cellfun(@(x)graph.prob(x), neighboursIdx, 'uni', 0);
graph.neighBorderIdx = cellfun(@(x)graph.borderIdx(x), neighboursIdx, 'uni', 0);
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
segmentMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta.mat');
segmentMeta.point = segmentMeta.point';
segmentMeta = connectEM.addSegmentClassInformation(p, segmentMeta);

%%
workingFolder = '/gaba/scratch/mberning/aggloGridSearch2/';
runs = dir([workingFolder 'search01_*']);
for i=1:length(runs)
    results{i} = dir([workingFolder runs(i).name filesep '*.mat']);
end

%%
load([workingFolder 'parameters.mat']);
load([workingFolder 'state.mat']);

endingsJoined = cellfun(@(x)regexp(x, '# endings \(all above\): (\d*)\n', 'tokens'), state, 'uni', 0);
endingsJoined = cellfun(@(x)cellfun(@str2double, x), endingsJoined, 'uni', 0);

figure; hold on;
for i=1:length(endingsJoined)
    plot(endingsJoined{i});
    label{i} = num2str(inputArgumentsAxons(i,:));
end
legend(label);


%% Result 76
% Closest to Kevin's parameters, 20 instead of 30 border threshold
find(all(bsxfun(@eq, inputArgumentsAxons, [0.7 0.9 0.8 20 0.5]),2))
load([workingFolder runs(76).name filesep results{76}(end-1).name], 'axonsNew');
y = connectEM.evaluateAggloMetaMeta(graph, axonsNew, [], 'dir_search01_run76', segmentMeta);

%% Result 41
% Closest to Kevin's parameters, 20 instead of 30 border threshold
find(all(bsxfun(@eq, inputArgumentsAxons, [0.8 0.95 0.7 40 0.3]),2))
load([workingFolder runs(44).name filesep results{44}(end-1).name], 'axonsNew');
y = connectEM.evaluateAggloMetaMeta(graph, axonsNew, [], 'dir_search01_run44', segmentMeta);
