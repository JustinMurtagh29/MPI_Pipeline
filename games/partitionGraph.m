function agglo = partitionGraph( graph, nrCluster )
    %CLUSTER adj matrix into nrCluster cluster

    % Remember all ids in edge array
    allId = unique(graph.edges(:));

    % Create adjaceny matrix with continous numbering
    edges = zeros(size(graph.edges));
    for i=1:length(allId)
        edges(graph.edges == allId(i)) = i;
    end
    adj = accumarray(edges, -log(graph.prob), ...
        [length(allId) length(allId)], @max, 0, 1);
    adj = adj + adj';

    % Find all shortest paths
    dist = graphallshortestpaths(adj);
    % Take only path length to a maximum of 1 (cut off long tail)
    dist(dist > 1) = 0;
    % Make double stochastic
    dist = dist ./ (1.0001*max(sum(dist)));  % make largest sum a little smaller
    % than 1 to make sure no entry of C becomes negative
    dist = dist + sparse(1:size(dist,1), 1:size(dist,2), 1-sum(dist));

    % Make symetric again as graphallshortestpaths yield A - A' ~ 10^-15
    dist = (dist + dist') ./ 2;

    % Check whether symmetric (due to generation above it should be)
    if ~issymmetric(dist)
        error('Matrix not symmetric');
    end

    % Check whether real (due to generation above it should be)
    if ~isreal(dist)
        error('Matrix not real');
    end

    % Calculate nrCluster largest eigenvalues of matrix
    [V,D] = eig(dist, 'nobalance');
    chosenV = V(:,(end-nrCluster+1):end);
    chosenD = diag(D((end-nrCluster+1):end,(end-nrCluster+1):end));

    if any(chosenD < eps)
        warning('Very small eigenvalues used');
    end

    % K-means clustering
    clusterIdx = kmeans(chosenV, nrCluster, 'Distance', 'cosine', 'Start', 'sample', 'Replicates', 1);

    % Make everything nice and shiny for output
    [sortedIdx, perm] = sort(clusterIdx);
    allIdSorted = allId(perm);
    clusterLengths = diff([0; find(diff(sortedIdx) > 0); length(sortedIdx)]);
    agglo = mat2cell(allIdSorted, clusterLengths);

end

