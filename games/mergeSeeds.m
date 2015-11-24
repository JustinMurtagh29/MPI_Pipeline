function newSeeds = mergeSeeds(seeds)
    % Agglomerates supervoxel given the graph, seeds and a probability threshold

    display('Merging redunancies');
    allSeeds = cat(1,seeds{:}); 
    maxId = max(allSeeds);
    % Preallocate space 
    seedSize = cellfun(@length, seeds);
    nrElements = sum(seedSize.^2);
    adj = spalloc(maxId,maxId,nrElements);
    tic;
    for i=1:length(seeds)
        [p,q] = meshgrid(seeds{i}, seeds{i});
        adj(p(:),q(:)) = 1;
        Util.progressBar(i, length(seeds));
    end
    newSeeds = Graph.findConnectedComponents(adj, false, true);

end

