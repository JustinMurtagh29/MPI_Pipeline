
info = Util.runInfo(false);

load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');

if ~exist('heuristics','var')
    heuristics = load(fullfile(p.saveFolder, 'heuristicResult.mat'),'myelinScore');
end
if ~exist('graph','var')
    graph = load(fullfile(p.saveFolder, 'graphNew.mat'),'edges');
end

% get all edges between myelinated segments (+ selfEdges)
myelinThresh = 0.5;
myelinEdges = cat(1,repmat(find(heuristics.myelinScore>myelinThresh),1,2),graph.edges(all(heuristics.myelinScore(graph.edges)>myelinThresh,2),:));

% calc connected components
myelinAgglos = Graph.findConnectedComponents(myelinEdges,0,1);

save(fullfile(p.saveFolder,'aggloState','myelinAgglos.mat'),'myelinAgglos','info')