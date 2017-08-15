for idx_top = 1 : 20
    idx_top
    m = load(['/tmpscratch/kboerg/chiasmarun/resultchiasma_' num2str(ind(idx_top))]);
    p.sphereRadiusInner=4000;
    p.sphereRadiusOuter=Inf;
    [thisNodes, thisEdges, thisProb] = ...
        connectEM.detectChiasmataPruneToSphere( ...
        nodesScaled, temp{which_col}.output.edges, ...
        ones(size(temp{which_col}.output.edges,1),1), p, temp{which_col}.output.ccCenterIdx(m.i_here));
    mask= ~ismember(m.thisNodes,thisNodes,'rows'); 
    for idx =1 : length(m.C)
        m.C2{idx} =intersect(m.C{idx},find(mask));
    end
    connectEM.generateSkeletonFromAgglo(m.thisEdges, round(bsxfun(@times,m.thisNodes,1./[11.24,11.24,28])), m.C2, arrayfun(@(x) sprintf('skel_%d_%d',ind(idx_top),x),1:numel(m.C),'uni',0), './', max(skelsegids{1}));
    copyfile(ff.filenames{fflookup(idxtemp{which_col}(ind(idx_top)))}, '.')
end