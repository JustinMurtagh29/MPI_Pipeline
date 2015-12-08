function partition = splitPartitionByCC(graph, partition)

    problemsFound = 0;
    for i=1:length(partition)
        idx = all(ismember(graph.edges, partition{i}),2);
        edges = graph.edges(idx,:);
        cc = Graph.findConnectedComponents(edges, false, true);
        if length(cc) > 1
            partition{i} = [];
            for j=1:length(cc)
                partition{end+1} = cc{j};
            end
            problemsFound = problemsFound + 1;
        end
    end
    display(['Found ' num2str(problemsFound) ' components in partition which were not connected']);

end

