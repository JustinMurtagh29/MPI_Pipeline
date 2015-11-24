function newSeeds = agglomerateSG(graph, seeds, threshold, steps)
    % Agglomerates supervoxel given the graph, seeds and a probability threshold

    display('Generating efficent thresholded graph representation');
    tic;
    % Determine unique IDs in non-thresholded graph
    allId = unique(graph.edges(:));
    % Create adjaceny matrices (including self edges at allId)
    edges = graph.edges(graph.prob > threshold,:);
    adj = sparse([edges(:,1); edges(:,2); allId], [edges(:,2); edges(:,1); allId], 1);
    adjS = adj^steps;
    toc;

    display('Agglomerating supervoxel');
    newSeeds = cell(size(seeds));
    tic;
    for i=1:length(seeds)
        % Generate logical vector representation of seed
        thisSeed = sparse(size(adj,1),1);
        thisSeed(seeds{i}) = 1;
        % Multiply with matrix to power of steps (precalculated above)
        temp = adjS*thisSeed;
        newSeeds{i} = find(temp);
        Util.progressBar(i, length(seeds));
    end

end

