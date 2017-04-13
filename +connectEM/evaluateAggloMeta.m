function y = evalutateAggloMeta(skelpath, graph, segmentMeta, agglos, p)
y.recall_col = {};
y.splits_col = [];
y.mergers_col = [];
y.length_col = [];
y.covered_col = [];
y.numagglomerates_col = [];
for file_idx = 1 : length(skelpath)
    skel = skeleton(skelpath{file_idx});

    for idx = 1 : length(skel.nodes)
        skel.nodes{idx}(:, 1 : 3) = bsxfun(@minus, skel.nodes{idx}(:, 1 : 3),  [1195, 1515, 115]-(129 - [25 25 10]));

        % skel = connectEM.evaluateAggloCleanSkel(skel, idx, any(skel.nodes{idx} <= 0, 2));
        ids = Seg.Global.getSegIds(p, skel.nodes{idx}(:, 1 : 3));
        skel = connectEM.evaluateAggloCleanSkel(skel, idx, ids == 0);
        %skel = connectEM.intermediateNodes(skel, idx, 400)
        [recall, splits, mergers, validnodes, numagglomerates] = connectEM.evaluateAgglo(agglos, segmentMeta.point', skel, idx, ids, graph.neighbours);
        y.recall_col{end + 1} = recall;
        y.splits_col(end + 1) = splits;
        y.mergers_col(end + 1) = mergers;
        y.length_col(end + 1) = full(sum(sum(skel.createWeightedAdjacencyMatrix(idx)))) / 2;
        good_edges = skel.edges{idx}(all(ismember(skel.edges{idx}, validnodes), 2), :);
        y.covered_col(end + 1) = sum(sqrt(sum(bsxfun(@times, skel.nodes{idx}(good_edges(:, 1), 1 : 3) - skel.nodes{idx}(good_edges(:, 2), 1 : 3), [11.24, 11.24, 28]).^2, 2)));
        y.numagglomerates_col(end + 1) = numagglomerates;
    end
end
