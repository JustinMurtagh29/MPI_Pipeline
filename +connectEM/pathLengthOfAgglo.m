function pL = pathLengthOfAgglo(segmentMeta, agglo, voxelSize);
    if ~exist('voxelSize','var') || isempty(voxelSize)
        voxelSize = [11.24, 11.24, 28];
    end
    pos = segmentMeta.point(agglo,:);
    pos = bsxfun(@times, pos, voxelSize);
    adj = squareform(pdist(pos));
    adj(adj > 5000) = 0;
    tree = graphminspantree(sparse(adj), 'Method', 'Kruskal');
    [edges(:,1), edges(:,2)] = find(tree);
    edgeDiff = pos(edges(:,1),:)-pos(edges(:,2),:);
    pL = sqrt(sum(edgeDiff.^2,2));
    pL = sum(pL)/1000;

end

