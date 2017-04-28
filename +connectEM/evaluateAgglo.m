function [recall, splits, mergers, validnodes, foundAgglomerates, connM] = evaluateAgglo(agglomerates, segmentMeta, skel, skelAsIds, neighbours, limitaggloNum, limitaggloSize, agglos_reverse, maxTube)
    % Occurences of each agglomerate
    foundAgglomeratesPre = agglos_reverse(skelAsIds(skelAsIds ~= 0));
    foundAgglomeratesPre(foundAgglomeratesPre == 0) = [];
    % GT node per agglo threshold
    if ~isempty(foundAgglomeratesPre)
        foundAgglomerates = feval(@(x)x(x(:, 2) > limitaggloNum, 1), tabulate(foundAgglomeratesPre));
    else
        foundAgglomerates = [];
    end
    % Volume per agglo threshold
    aggloSize = cellfun(@(x)sum(segmentMeta.voxelCount(x)), agglomerates(foundAgglomerates));
    foundAgglomerates(aggloSize < limitaggloSize) = [];
    % Metrics
    recall = [length(intersect(cell2mat(agglomerates(foundAgglomerates)), skelAsIds(skelAsIds ~= 0))), length(unique(skelAsIds(skelAsIds ~= 0)))];
    validnodes = find(ismember(skelAsIds, cell2mat(agglomerates(foundAgglomerates))));
    % Merger calculation
    mergers = 0;
    scalize = @(x)bsxfun(@times,x,[11.24, 11.24, 28]);
    for idx = 1 : length(foundAgglomerates)
        points = segmentMeta.point(agglomerates{foundAgglomerates(idx)}, :);
        if size(points, 1) > 1000
            points = points(randperm(size(points, 1), 1000), :);
        end
        if max(min(pdist2(scalize(points), scalize(segmentMeta.point(skelAsIds(skelAsIds > 0), :))), [], 2)) > maxTube
            mergers = mergers + 1;
        end
    end
    % Connectivity matrix of found agglomerates (which touch ...)
    connM = zeros(length(foundAgglomerates));
    %{
    for idx1 = 1 : length(foundAgglomerates)
        for idx2 = idx1 : length(foundAgglomerates)
            if any(cellfun(@(x)any(ismember(x, agglomerates{foundAgglomerates(idx2)})), neighbours(agglomerates{foundAgglomerates(idx1)})))
                connM(idx2, idx1) = 1; %only set the lower diagonal matrix
            end
        end
    end
    %}
    % Split calculation: Not currently in use
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
