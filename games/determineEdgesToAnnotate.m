function [edgesToAnnotate, probEdges, edgesToForget, probEdges2] = determineEdgesToAnnotate(graph, unusedIds)
    % Determine which edges to annotate in current cube
    
    % First all edges to unused IDs with highest probability
    for i=1:length(unusedIds)
        idxE = any(ismember(graph.edges, unusedIds(i)),2);
        [mp,idxMP] = max(graph.prob(idxE));
        idxE = find(idxE);
        edgesToAnnotate(i,:) = graph.edges(idxE(idxMP),:);
        probEdges(i) = mp;
    end
    % Split by percentage
    idxF = probEdges < .4;
    edgesToForget = edgesToAnnotate(idxF,:);
    probEdges2 = probEdges(idxF);
    edgesToAnnotate(idxF,:) = [];
    probEdges(idxF) = [];

end

