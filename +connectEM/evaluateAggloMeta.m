function y = evalutateAggloMeta(skelpath, graph, segmentMeta, agglos, p, nameOfRun, mainFolder, limitaggloNum, limitaggloSize)
segmentMeta2 = segmentMeta;
segmentMeta2.point = segmentMeta2.point'
doc_folder = [mainFolder, nameOfRun, '/'];
mkdir(doc_folder);
for file_idx = 1 : length(skelpath)
    skel = skeleton(skelpath{file_idx});
    assert(length(skel.nodes) == 1)
    idx = 1;
    if isempty(strfind(skel.parameters.experiment.name, 'ROI2017'))
        skel.nodes{idx}(:, 1 : 3) = bsxfun(@minus, skel.nodes{idx}(:, 1 : 3),  [1195, 1515, 115]-(129 - [25 25 10]));
    end
    skel.nodesNumDataAll{idx}(:, 3 : 5)= skel.nodes{idx}(:, 1 : 3);
    % skel = connectEM.evaluateAggloCleanSkel(skel, idx, any(skel.nodes{idx} <= 0, 2));
    % skel = connectEM.intermediateNodes(skel, idx, 400);
    assert(graphconncomp(skel.createAdjacencyMatrix(idx)) == 1);
    assert(~skel.hasCircle(idx));
    ids = Seg.Global.getSegIds(p, skel.nodes{idx}(:, 1 : 3));
    skel = connectEM.evaluateAggloCleanSkel(skel, idx, ids == 0);
    ids = Seg.Global.getSegIds(p, skel.nodes{idx}(:, 1 : 3));
    [recall, splits, mergers, validnodes, foundAgglomerates, connM] = connectEM.evaluateAgglo(agglos, segmentMeta2, skel, idx, ids, graph.neighbours, limitaggloNum, limitaggloSize);
    y.recall_col{file_idx} = recall;
    y.splits_col(file_idx) = splits;
    y.mergers_col(file_idx) = mergers;
    y.length_col(file_idx) = full(sum(sum(skel.createWeightedAdjacencyMatrix(idx)))) / 2;
    good_edges = skel.edges{idx}(all(ismember(skel.edges{idx}, validnodes), 2), :);
    y.covered_col(file_idx) = sum(sqrt(sum(bsxfun(@times, skel.nodes{idx}(good_edges(:, 1), 1 : 3) - skel.nodes{idx}(good_edges(:, 2), 1 : 3), [11.24, 11.24, 28]).^2, 2)));
    y.foundAgglomerates_col{file_idx} = foundAgglomerates;

    y.connM{file_idx} = connM;
    connectEM.skeletonFromAgglo(graph.edges, segmentMeta2, agglos(foundAgglomerates), num2str(file_idx), doc_folder)
    skel.write([doc_folder, 'skel_' num2str(file_idx)  '.nml']);
end
