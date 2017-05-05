folder = '/gaba/u/kboerg/code/pipeline/+connectEM/evaluationData/new_axon_gt/'
liste = dir([folder '*.nml']);
for idx = 1 : 10
    skel = skeleton([folder liste(idx).name]);
    agglos = gridAgglo_05{564}.metrics.axon1.foundAgglomerates_col{idx};
    store(idx) = 0;
    for idx2 = 1 : length(agglos)
        [Tree, pred] = graphminspantree(sparse(squareform(pdist(bsxfun(@times, segmentMeta.point(gridAgglo_05{564}.axonsFinal{agglos(idx2)}, :), [11.24, 11.24, 28])))));
        store(idx) = store(idx) + sum(Tree(:));
    end
    store(idx) = store(idx) / (full(sum(sum(skel.createWeightedAdjacencyMatrix(1))))/2) /gridAgglo_05{564}.metrics.axon1.recall_col{idx}(1)*gridAgglo_05{564}.metrics.axon1.recall_col{idx}(2);
end
