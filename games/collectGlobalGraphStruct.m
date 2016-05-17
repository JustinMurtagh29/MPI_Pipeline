function collectGlobalGraphStruct (p)

load([p.saveFolder 'globalEdges.mat']);
load([p.saveFolder 'globalCoMList.mat']);
load([p.saveFolder 'globalGPProbList.mat']);

[ neighbours, neighboursIdx ] = Graph.edges2Neighbors(edges);

neighProb = cellfun(@(x) prob(x),neighboursIdx,'UniformOutput',false);

graph.neighbours = neighbours;
graph.neighProb = neighProb;
graph.edges = edges;
graph.prob = prob;

save([p.saveFolder 'graph.mat'],'graph','-v7.3');

end
