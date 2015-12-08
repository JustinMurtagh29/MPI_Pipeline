function agglo = partitionGraph( graph, nrCluster )
    %Cluster adj matrix into nrCluster cluster

    % Remeber all ids in edge array
    allId = unique(graphR.edges(:));

    % Create adjaceny matrix with continous numbering
    edges = zeros(size(graphR.edges));
    for i=1:length(allId)
        edges(graphR.edges == allId(i)) = i;
    end
    adj = accumarray(edges, graphR.prob, ...
        [length(allId) length(allId)], @max);
    adj = adj + adj';

    % Normalize matrix
    adj = adj ./ (1.001*max(sum(adj)));
    adj = adj + diag(1-sum(adj));

    % Calculate nrCluster largest eigenvalues of matrix
    [eigenVectors, eigenValues] = eig(adj, 'nobalance', 'vector');
    [~,perm] = sort(eigenValues, 'descend');
    eigenVectors = eigenVectors(:,perm(1:nrCluster));
    eigenValues = eigenValues(perm(1:nrCluster));

    % K-means clustering
    clusterIdx = kmeans(eigenVectors, nrCluster, ...
        'Distance', 'cosine', 'Start', 'sample', 'Replicates', 1);

    % Make everything nice and shiny for output
    [sortedIdx, perm] = sort(clusterIdx);
    allIdSorted = allId(perm);
    clusterLengths = diff([0; find(diff(sortedIdx) > 0); length(sortedIdx)]);
    agglo = mat2cell(allIdSorted, clusterLengths);

end

