function y = evalutateAggloMeta(skelpath, graph, segmentMeta, agglos, p)
y.recall_col = {};
y.splits_col = [];
y.mergers_col = [];
for file_idx = 1 : length(skelpath)
    skel = skeleton(skelpath{file_idx});

    for idx = 1 : length(skel.nodes)
        skel.nodes{idx}(:, 1 : 3) = bsxfun(@minus, skel.nodes{idx}(:, 1 : 3),  [1195, 1515, 115]-(129 - [25 25 10]));

        % skel = connectEM.evaluateAggloCleanSkel(skel, idx, any(skel.nodes{idx} <= 0, 2));
        ids = Seg.Global.getSegIds(p, skel.nodes{idx}(:, 1 : 3));
        skel = connectEM.evaluateAggloCleanSkel(skel, idx, ids == 0);
        %skel = connectEM.intermediateNodes(skel, idx, 400)
        [recall, splits, mergers] = connectEM.evaluateAgglo(agglos, segmentMeta.point', skel, idx, ids, graph.neighbours);
        y.recall_col{end + 1} = recall;
        y.splits_col(end + 1) = splits;
        y.mergers_col(end + 1) = mergers;
    end
end
