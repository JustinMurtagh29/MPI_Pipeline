function [recall, splits, mergers, validnodes, foundAgglomerates, connM] = evaluateAgglo(agglomerates, segmentMeta2, skel, skelIdx, skelAsIds, neighbours, limitaggloNum, limitaggloSize, agglos_reverse)
    maxTube = 10000;
    foundAgglomeratesPre = setdiff(agglos_reverse(intersect(skelAsIds(skelAsIds ~= 0), 1:agglos_reverse)), 0);
    if ~isempty(foundAgglomeratesPre)
        foundAgglomerates = feval(@(x)x(x(:, 2) > limitaggloNum, 1), tabulate(foundAgglomeratesPre));
    else
        foundAgglomerates = [];
    end
    aggloSize = cellfun(@(x)sum(segmentMeta2.point(x)), agglomerates(foundAgglomerates));
    foundAgglomerates(aggloSize < limitaggloSize) = [];
    recall = [length(intersect(cell2mat(agglomerates(foundAgglomerates)), skelAsIds(skelAsIds ~= 0))), length(unique(skelAsIds(skelAsIds ~= 0)))];
    validnodes = find(ismember(skelAsIds, cell2mat(agglomerates(foundAgglomerates))));
    mergers = 0;
    scalize = @(x)bsxfun(@times,x,[11.24, 11.24, 28]);
    for idx = 1 : length(foundAgglomerates)
        if max(min(pdist2(scalize(segmentMeta2.point(agglomerates{foundAgglomerates(idx)}, :)), scalize(segmentMeta2.point(skelAsIds(skelAsIds > 0), :))), [], 2)) > maxTube
            mergers = mergers + 1;
        end
    end
    connM = zeros(length(foundAgglomerates));
    for idx1 = 1 : length(foundAgglomerates)
        for idx2 = idx1 : length(foundAgglomerates)
            if any(cellfun(@(x)any(ismember(x, agglomerates{foundAgglomerates(idx2)})), neighbours(agglomerates{foundAgglomerates(idx1)})))
                connM(idx2, idx1) = 1; %only set the lower diagonal matrix
            end
        end
    end
    if graphconncomp(sparse(connM), 'Directed', false) ~= 1
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
            distances = graphallshortestpaths(sparse(connM), 'Directed', false);
            [~, idx2] = min(distances(inThere, endpointsMeta(idx)));
            [~, path] = graphshortestpath(sparse(connM), endpointsMeta(idx), endpointsMeta(idx2));
            inThere = [inThere, path];
            splits = splits + length(path);
        end
    end
end
