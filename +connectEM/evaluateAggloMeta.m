function y = evalutateAggloMeta(skelpath, graph, segmentMeta, agglos, p, nameOfRun, mainFolder)
segmentMeta2 = segmentMeta;
segmentMeta2.point = segmentMeta2.point'
doc_folder = [mainFolder, nameOfRun, '/'];
mkdir(doc_folder);
for file_idx = 1 : length(skelpath)
    skel = skeleton(skelpath{file_idx});
    assert(length(skel.nodes) == 1)
    idx = 1;
    skel.nodes{idx}(:, 1 : 3) = bsxfun(@minus, skel.nodes{idx}(:, 1 : 3),  [1195, 1515, 115]-(129 - [25 25 10]));
    skel.nodesNumDataAll{idx}(:, 3 : 5)= skel.nodes{idx}(:, 1 : 3);
    % skel = connectEM.evaluateAggloCleanSkel(skel, idx, any(skel.nodes{idx} <= 0, 2));
    ids = Seg.Global.getSegIds(p, skel.nodes{idx}(:, 1 : 3));
    skel = connectEM.evaluateAggloCleanSkel(skel, idx, ids == 0);
    %skel = connectEM.intermediateNodes(skel, idx, 400)
    [recall, splits, mergers, validnodes, foundAgglomerates] = connectEM.evaluateAgglo(agglos, segmentMeta.point', skel, idx, ids, graph.neighbours, 2);
    y.recall_col{file_idx} = recall;
    y.splits_col(file_idx) = splits;
    y.mergers_col(file_idx) = mergers;
    y.length_col(file_idx) = full(sum(sum(skel.createWeightedAdjacencyMatrix(idx)))) / 2;
    good_edges = skel.edges{idx}(all(ismember(skel.edges{idx}, validnodes), 2), :);
    y.covered_col(file_idx) = sum(sqrt(sum(bsxfun(@times, skel.nodes{idx}(good_edges(:, 1), 1 : 3) - skel.nodes{idx}(good_edges(:, 2), 1 : 3), [11.24, 11.24, 28]).^2, 2)));
    y.numagglomerates_col(file_idx) = length(foundAgglomerates);
    connectEM.skeletonFromAgglo(graph.edges, segmentMeta2, agglos(foundAgglomerates), num2str(file_idx), doc_folder)
    skel.write([doc_folder, 'skel_' num2str(file_idx)  '.nml']);
end
