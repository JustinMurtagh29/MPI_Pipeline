function [graphR, unusedIds] = restrictSGtoAboveProbability(graph, lowerT, upperT)

    edgeIdx = graph.prob > lowerT & graph.prob < upperT;
    graphR.edges = graph.edges(edgeIdx,:);
    graphR.prob = graph.prob(edgeIdx);
    graphR.cubeLI = graph.cubeLI(edgeIdx);
    graphR.borderCentroid = graph.borderCentroid(edgeIdx,:);

    allIdBefore = unique(graph.edges);
    allIdAfter = unique(graphR.edges);
    unusedIds = setdiff(allIdBefore, allIdAfter);

end

