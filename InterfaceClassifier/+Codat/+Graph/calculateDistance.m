function d = calculateDistance( graph, startNode, maxDist )
%CALCULATEDISTANCE Calculate all nodes within a maximal distance
%from a starting node on an undirected graph.
% INPUT graph: Either a sparse symmetric adjacency matrix or one direction
%              edge list.
%       startNode: ID of the start node.
%       maxDist: Maximal distance to consider in number of edges.
% OUTPUT d: Vector of length number of nodes specifying the distance for
%           each node from the target node. If the distance exceeds maxDist
%           then -1 is returned.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

% Create sparse undirected adjaceny matrix if not supplied
if ~issparse(graph)
    graph = double(graph);
    maxValue = max(graph(:));
    graph = sparse(graph(:,1),graph(:,2),1,maxValue,maxValue);
end

%make graph symmetric
graph = graph + graph';

d = -ones(size(graph,1),1);
d(startNode) = 0;
nodes = false(size(graph,1),1);
nodes(startNode) = 1;
nodesVisited = nodes;

for i = 1:maxDist
    nodes = graph*nodes;
    nodes(nodesVisited) = false;
    nodesVisited = nodesVisited | nodes;
    d(d < 0 & nodes) = i;
end


end

