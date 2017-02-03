function partition = partitionWholeDataset(graph, threshold)
    % Do a simple partition by thresholding graph and calculating CC
    
    remainingEdges = graph.edges(graph.prob > threshold,:);
    partition = Graph.findConnectedComponents(remainingEdges, true, true);  

end


