function pL = pathLengthOfAgglo(segmentMeta, agglo);

    pos = segmentMeta.point(agglo,:);
    pos = bsxfun(@times, pos, [11.24 11.24 28]);
    adj = squareform(pdist(pos));
    adj(adj > 5000) = 0;
    tree = graphminspantree(sparse(adj), 'Method', 'Kruskal');
    [edges(:,1), edges(:,2)] = find(tree);
    edgeDiff = pos(edges(:,1),:)-pos(edges(:,2),:);
    pL = sqrt(sum(edgeDiff.^2,2));
    pL = sum(pL)/1000;

end

