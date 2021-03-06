function getAxonPathLength(ID,param,superAgglos,batchBoundaries,dataDir)

superAgglos = superAgglos(batchBoundaries(ID):batchBoundaries(ID+1)-1);

for i=1:length(superAgglos)
    nodes = superAgglos(i).nodes(:,1:3);
    edges = minimalSpanningTree(nodes);
    distances = [];
    for j=1:size(edges,1)
        distances(j,1) = pdist(bsxfun(@times,[nodes(edges(j,1),:);nodes(edges(j,2),:)], [11.24 11.24 28]));
    end
    aggloLength(i,1) = sum(distances);
end

directory = strcat(dataDir,'batch',int2str(ID),'.mat');
save(directory, 'aggloLength')
end

function edges = minimalSpanningTree(com)
    if size(com,1) < 2
        edges = [];
    else
        % Minimal spanning tree
        adj = squareform(pdist(bsxfun(@times, com, [11.24 11.24 28])));
        tree = graphminspantree(sparse(adj), 'Method', 'Kruskal');
        [edges(:,1), edges(:,2)] = find(tree);
    end
end

