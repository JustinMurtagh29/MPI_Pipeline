function [splits, mergers] = evaluate(agglomerates, CoMs, skel, skelIdx, skelAsIds, graphM)
    maxTube = 10000
    allAglomerates = find(ismember(agglomerates, skelsAsIds));
    mergers = 0;
    scalize = @(x)bsxfun(@times,x,[11.24, 11.24, 28]);
    for idx = 1 : length(allAglomerates)
        if min(max(pdist2(scalize(CoMs(allAglomerates, :)), scalize(CoMs(skelAsIds, :))))) > maxTube
            mergers = mergers + 1;
        end
    end
    metaConnectivityMatrix = zeros(length(allAglomerates));
    for idx1 = 1 : length(allAglomerates)
        for idx2 = idx1 : length(allAglomerates)
            assert(graphM(1, 1)) % assumes that the diagonal is set in graphM
            if any(graphM(allAglomerates{idx1}, allAglomerates{idx2}))
                metaConnectivityMatrix(idx2, idx1) = 1; %only set the lower diagonal matrix
            end
        end
    end
    if graphconncomp(metaConnectivtyMatrix, 'Directed', false) ~= 1
        disp('skeleton does not stick together')  % todo interpolate nodes
        splits = -1;
        return;
    end
    endpoints = skelAsIds(sum(skel.createAdjacencyMatrix(1))==1);
    endpointsMeta = find(ismember(endpoints, agglomerates));
    inThere = [endpointsMeta(1)];
    splits = 0;
    for idx = 2 : length(endpointsMeta)
        distances = graphallshortestpath(metaConnectivityMatrix, 'Directed', false)
        [~, idx2] = min(distances(inThere, endpointsMeta(idx)));
        [~, path] = graphshortestpath(metaConnectivityMatrix, endPointsMeta(idx), endPointsMeta(idx2));
        inThere = [inThere, path]
        splits = splits + length(path);
    end
end
