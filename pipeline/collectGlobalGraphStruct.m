function collectGlobalGraphStruct (p)

load([p.saveFolder 'globalEdges.mat']);
load([p.saveFolder 'globalCoMList.mat']);
load([p.saveFolder 'globalGPProbList.mat']);
load([p.saveFolder 'globalBorder.mat']);

% Correspondences were missing, load here for now, decide whether to also save globally
correspondences = Seg.Global.getGlobalCorrespondences(p);
edges = cat(1, edges, correspondences);
prob = cat(1, prob, ones(length(correspondences),1));

[ neighbours, neighboursIdx ] = Graph.edges2Neighbors(edges);

neighProb = cellfun(@(x) prob(x),neighboursIdx,'UniformOutput',false);

graph = struct;
graph.neighbours = neighbours;
graph.neighProb = neighProb;
graph.edges = edges;
graph.prob = prob;
graph.borderCentroid = borderCoM;
graph.borderSize = borderSize;
graph.borderArea = borderArea;

save([p.saveFolder 'graph.mat'],'-struct','graph','-v7.3');

end





