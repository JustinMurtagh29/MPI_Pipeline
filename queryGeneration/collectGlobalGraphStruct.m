%{
function collectGlobalGraphStruct (p)

load([p.saveFolder 'globalEdges.mat']);
load([p.saveFolder 'globalCoMList.mat']);
load([p.saveFolder 'globalGPProbList.mat']);
load([p.saveFolder 'globalBorder.mat']);


[ neighbours, neighboursIdx ] = Graph.edges2Neighbors(edges);

neighProb = cellfun(@(x) prob(x),neighboursIdx,'UniformOutput',false);

graph.neighbours = neighbours;
graph.neighProb = neighProb;
graph.edges = edges;
graph.prob = prob;
graph.borderCentroid = borderCoM;
graph.borderSize = borderSize;
graph.borderArea = borderArea;

save([p.saveFolder 'graph.mat'],'graph','-v7.3');

end
%}


function collectGlobalGraphStruct (p)

load([p.saveFolder 'globalEdges.mat']);
load([p.saveFolder 'globalCoMList.mat']);
load([p.saveFolder 'globalGPProbList.mat']);
load([p.saveFolder 'globalBorder.mat']);


[ neighbours, neighboursIdx ] = Graph.edges2Neighbors(edges);

neighProb = cellfun(@(x) prob(x),neighboursIdx,'UniformOutput',false);

graph.neighbours = neighbours;
graph.neighProb = neighProb;
graph.edges = edges;
graph.prob = prob;
graph.borderCentroid = borderCoM;
graph.borderSize = borderSize;
graph.borderArea = borderArea;

save([p.saveFolder 'graph.mat'],'graph','-v7.3');

end





