function [partition, remainingEdges] = partitionWholeDataset(graph, threshold)
    % Do a simple partition by thresholding graph and calculating CC
    
    remainingEdges = graph.edges(graph.prob > threshold,:);
    [partition,~,remainingEdges] = Graph.findConnectedComponents(remainingEdges, true, true);  

end


