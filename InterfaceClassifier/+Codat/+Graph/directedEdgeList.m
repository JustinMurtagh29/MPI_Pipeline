function edges = directedEdgeList( edges, rootNode )
%DIRECTEDEDGELIST Create a directed graph by starting from a root node.
% INPUT edges: [Nx2] integer array. Each row defines an edge between the
%              corresponding nodes. (The entries in edges can be indices of
%              nodes in a skel or IDs of nodes).
%       rootNode: Node ID in edges which to set as the root node.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

nodeList = rootNode;
isDirected = edges(:,1) == rootNode;

%start from root node and always directe eges which have not yet been
%directed away from the current node
while ~isempty(nodeList)
    currNode = nodeList(1);
    wrongDir= edges(:,2) == currNode & ~isDirected;
    edges(wrongDir,[2, 1]) = edges(wrongDir,:);
    isDirected(wrongDir) = true;
    isDirected(edges(:,1) == currNode) = true;
    nodeList = cat(1,nodeList,edges(edges(:,1) == currNode,2));
    nodeList(1) = [];
end

end

