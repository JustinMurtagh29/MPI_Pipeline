function [recall, splits, mergers] = evaluateAgglo(agglomerates, CoMs, skel, skelIdx, skelAsIds, neighbours)
    maxTube = 10000;
    allAglomerates = find(cellfun(@(x)any(ismember(x, skelAsIds(skelAsIds ~= 0))), agglomerates));
    recall = [length(intersect(cell2mat(agglomerates(allAglomerates)), skelAsIds(skelAsIds ~= 0))), length(skelAsIds(skelAsIds ~= 0))]
    mergers = 0;
    scalize = @(x)bsxfun(@times,x,[11.24, 11.24, 28]);
    for idx = 1 : length(allAglomerates)
        if min(max(pdist2(scalize(CoMs(aglomerates{allAglomerates(idx)}, :)), scalize(CoMs(skelAsIds(skelAsIds > 0), :))))) > maxTube
            mergers = mergers + 1;
        end
    end
    metaConnectivityMatrix = zeros(length(allAglomerates));
    for idx1 = 1 : length(allAglomerates)
        for idx2 = idx1 : length(allAglomerates)
            if any(cellfun(@(x)any(ismember(x, agglomerates{allAglomerates(idx2)})), neighbours(agglomerates{allAglomerates(idx1)})))
                metaConnectivityMatrix(idx2, idx1) = 1; %only set the lower diagonal matrix
            end
        end
    end
    if graphconncomp(sparse(metaConnectivityMatrix), 'Directed', false) ~= 1
        disp('skeleton does not stick together')  % todo interpolate nodes
        splits = -1;
        return;
    end
    endpoints = skelAsIds(sum(skel.createAdjacencyMatrix(1))==1); % we assured in the outer function that endpoints don't have seg id 0
    endpointsMeta = find(cellfun(@(x)any(ismember(x, endpoints)), agglomerates));
    inThere = [endpointsMeta(1)];
    splits = 0;
    for idx = 2 : length(endpointsMeta)
        distances = graphallshortestpath(sparse(metaConnectivityMatrix), 'Directed', false)
        [~, idx2] = min(distances(inThere, endpointsMeta(idx)));
        [~, path] = graphshortestpath(sparse(metaConnectivityMatrix), endPointsMeta(idx), endPointsMeta(idx2));
        inThere = [inThere, path]
        splits = splits + length(path);
    end
end
