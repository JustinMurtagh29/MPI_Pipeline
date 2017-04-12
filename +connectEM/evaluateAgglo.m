function [recall, splits, mergers] = evaluateAgglo(agglomerates, CoMs, skel, skelIdx, skelAsIds, neighbours)
    maxTube = 10000;
    allAgglomerates = find(cellfun(@(x)any(ismember(x, skelAsIds(skelAsIds ~= 0))), agglomerates));
    recall = [length(intersect(cell2mat(agglomerates(allAgglomerates)), skelAsIds(skelAsIds ~= 0))), length(skelAsIds(skelAsIds ~= 0))];
    mergers = 0;
    scalize = @(x)bsxfun(@times,x,[11.24, 11.24, 28]);
    for idx = 1 : length(allAgglomerates)
        if min(max(pdist2(scalize(CoMs(agglomerates{allAgglomerates(idx)}, :)), scalize(CoMs(skelAsIds(skelAsIds > 0), :))))) > maxTube
            mergers = mergers + 1;
        end
    end
    metaConnectivityMatrix = zeros(length(allAgglomerates));
    for idx1 = 1 : length(allAgglomerates)
        for idx2 = idx1 : length(allAgglomerates)
            if any(cellfun(@(x)any(ismember(x, agglomerates{allAgglomerates(idx2)})), neighbours(agglomerates{allAgglomerates(idx1)})))
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
    endpointsMeta = find(cellfun(@(x)any(ismember(x, endpoints)), agglomerates(allAgglomerates)));
    if isempty(endpointsMeta)
        splits = -2;
    else
        inThere = [endpointsMeta(1)];
        splits = 0;
        for idx = 2 : length(endpointsMeta)
            distances = graphallshortestpaths(sparse(metaConnectivityMatrix), 'Directed', false);
            [~, idx2] = min(distances(inThere, endpointsMeta(idx)));
            [~, path] = graphshortestpath(sparse(metaConnectivityMatrix), endpointsMeta(idx), endpointsMeta(idx2));
            inThere = [inThere, path];
            splits = splits + length(path);
        end
    end
end
