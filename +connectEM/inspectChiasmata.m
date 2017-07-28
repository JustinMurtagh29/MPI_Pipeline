for idx_top = 1 : 5
    m = load(['/tmpscratch/kboerg/chiasmarun/resultchiasma_' num2str(ind(idx_top))]);
    mask= pdist2(m.thisNodes,nodesScaled(temp{which_col}.output.ccCenterIdx(m.i_here),:)) < 4000;
    mask2 = zeros(size(mask));
    mask2(mask) = 1 : sum(mask);
    for idx =1 : length(m.C)
        m.C{idx} = setdiff(mask2(m.C{idx}),0);
    end
    m.thisEdges=mask2(m.thisEdges);
    m.thisEdges(any(m.thisEdges==0,2),:) = [];
    m.thisNodes=m.thisNodes(mask,:);
     connectEM.generateSkeletonFromAgglo(m.thisEdges, round(bsxfun(@times,m.thisNodes,1./[11.24,11.24,28])), m.C, arrayfun(@(x) sprintf('skel_%d_%d',ind(idx_top),x),1:numel(m.C),'uni',0), './', max(skelsegids{1}));
    copyfile(ff.filenames{fflookup(idxtemp{which_col}(ind(idx_top)))}, '.')
end