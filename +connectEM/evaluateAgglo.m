function [recall, splits, mergers, validnodes, foundAgglomerates] = evaluateAgglo(agglomerates, CoMs, skel, skelIdx, skelAsIds, neighbours, limitagglo)
    maxTube = 10000;
    foundAgglomerates = find(cellfun(@(x)sum(ismember(skelAsIds(skelAsIds ~= 0), x)) > limitagglo, agglomerates));
    recall = [length(intersect(cell2mat(agglomerates(foundAgglomerates)), skelAsIds(skelAsIds ~= 0))), length(skelAsIds(skelAsIds ~= 0))];
    validnodes = find(ismember(skelAsIds, cell2mat(agglomerates(foundAgglomerates))));
    mergers = 0;
    scalize = @(x)bsxfun(@times,x,[11.24, 11.24, 28]);
    for idx = 1 : length(foundAgglomerates)
        if max(min(pdist2(scalize(CoMs(agglomerates{foundAgglomerates(idx)}, :)), scalize(CoMs(skelAsIds(skelAsIds > 0), :))), [], 2)) > maxTube
            mergers = mergers + 1;
        end
    end
    metaConnectivityMatrix = zeros(length(foundAgglomerates));
    for idx1 = 1 : length(foundAgglomerates)
        for idx2 = idx1 : length(foundAgglomerates)
            if any(cellfun(@(x)any(ismember(x, agglomerates{foundAgglomerates(idx2)})), neighbours(agglomerates{foundAgglomerates(idx1)})))
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
    endpointsMeta = find(cellfun(@(x)any(ismember(x, endpoints)), agglomerates(foundAgglomerates)));
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
