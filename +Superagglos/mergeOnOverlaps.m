function agglo = mergeOnOverlaps(aggloA, aggloB, varargin)
    % agglo = mergerOnOverlaps(aggloA, aggloB)
    %   This function takes a pair of super-agglomerates (`aggloA` and
    %   `aggloB`) and merges `aggloB` into `aggloA`. To avoid redundancies,
    %   `aggloB` is first stripped from all parts it has in common with
    %   `aggloA`.
    %
    % NOTE
    %   This function is not symmetric with respect to its input arguments.
    %   That is, `mergeOnOverlaps(A, B)` is not `mergeOnOverlaps(B, A)`.
    %   Typically, you'd want the "stronger" agglomerate to be `aggloA`.
    %
    % Options
    %   * scale: 1x3 vector which is used to scale the nodes before
    %     calculating the distance between the super-agglomerates.
    %     Default value: [1, 1, 1]
    %
    %   * overlapDistNm: scalar which determines the distance threshold (in
    %     nm) below which nodes from super-agglomerates A and B are
    %     considered to be overlapping.
    %     Default value: 100
    %
    %   * minLenNm: scalar which determines the length threshold (in nm)
    %     below which agglomerate stretches may be discarded.
    %     Default value: 2000
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    opts = struct;
    opts.scale = [1, 1, 1];
    opts.overlapDistNm = 100;
    opts.minLenNm = 2000;
    opts.use2015b = false;
    
    % apply user-specified values
    opts = Util.modifyStruct(opts, varargin{:});
    
    % sanity checks
    if opts.use2015b
        assert(all(aggloB.edges(:,2) >= aggloB.edges(:,1)));
    else
        assert(issorted(aggloB.edges, 2));
    end
    
    %% find and discard common nodes
    % calculate pair-wise distance
    distMat = pdist2( ...
        bsxfun(@times, aggloA.nodes(:, 1:3), opts.scale), ...
        bsxfun(@times, aggloB.nodes(:, 1:3), opts.scale));
    distVec = reshape(min(distMat, [], 1), [], 1);

    % discard common nodes
    nodeIdsB = find(distVec > opts.overlapDistNm);

    % grow out one step
    nodeIdsB = any(ismember(aggloB.edges, nodeIdsB), 2);
    nodeIdsB = unique(aggloB.edges(nodeIdsB, :));

    % discard edges with discarded nodes
   [~, relEdgesB] = ismember(aggloB.edges, nodeIdsB);
    relEdgesB = relEdgesB(all(relEdgesB, 2), :);
    
    %%
    % split into connected components
    adjMatB = sparse( ...
        relEdgesB(:, 2), relEdgesB(:, 1), 1, ...
        numel(nodeIdsB), numel(nodeIdsB));
   [numCompsB, lutB] = graphconncomp( ...
        adjMatB, 'Directed', false);
    
    % For each connected component of B we determine the edge which
    % connects it to super-agglomerate A by searching for the point of
    % closest proximity. This is only done if the B component has at least
    % the required length. Otherwise, the connecting edge is `[nan, nan]`.
    edgesAB = nan(numCompsB, 2);
    for curIdx = 1:numCompsB
        curNodeIds = nodeIdsB(lutB == curIdx);
        curLen = pathLen(bsxfun(@times, aggloB.nodes(curNodeIds, 1:3), opts.scale));

        % tiny component → ignore
        if curLen < opts.minLenNm; continue; end

        curDistMat = distMat(:, curNodeIds);
       [~, curMinIdx] = min(curDistMat(:));
       [curNodeA, curNodeB] = ind2sub(size(curDistMat), curMinIdx);

        edgesAB(curIdx, 1) = curNodeA;
        edgesAB(curIdx, 2) = curNodeIds(curNodeB);
    end

    compMaskB = ~isnan(edgesAB(:, 1));
    nodeIdsB = nodeIdsB(compMaskB(lutB));

    edgesAB = edgesAB(compMaskB, :);
   [~, edgesAB(:, 2)] = ismember(edgesAB(:, 2), nodeIdsB);
   [~, edgesB] = ismember(aggloB.edges, nodeIdsB);
    edgesB = edgesB(all(edgesB, 2), :);
    
    % globalize
    nodeCountA = size(aggloA.nodes, 1);
    edgesAB(:, 2) = edgesAB(:, 2) + nodeCountA;
    edgesB = edgesB + nodeCountA;
    
    %% building output agglomerate
    agglo = struct;
    agglo.nodes = cat(1, aggloA.nodes, aggloB.nodes(nodeIdsB, :));
    agglo.edges = cat(1, aggloA.edges, edgesAB, edgesB);
    
    if isfield(aggloA, 'endings') ...
            || isfield(aggloB, 'endings')
        agglo.endings = vertcat( ...
            aggloA.endings, nodeCountA + aggloB.endings);
        agglo.endings = setdiff(agglo.endings, edgesAB);
        agglo.endings = reshape(agglo.endings, [], 1);
    end
    
    if isfield(aggloA, 'solvedChiasma') ...
            || isfield(aggloB, 'solvedChiasma')
        agglo.solvedChiasma = vertcat( ...
            aggloA.solvedChiasma, aggloB.solvedChiasma);
        agglo.solvedChiasma(edgesAB) = false;
    end
    
    % sanity checks
    if opts.use2015b
        assert(all(agglo.edges(:,2) >= agglo.edges(:,1)));
    else
        assert(issorted(agglo.edges, 2));
    end
    
end

function len = pathLen(nodesNm)
    if size(nodesNm, 1) < 2
        len = 0;
        return;
    end
    
    adj = squareform(pdist(nodesNm));
    mst = graphminspantree(sparse(adj), 'Method', 'Kruskal');
    len = full(sum(mst(:)));
end