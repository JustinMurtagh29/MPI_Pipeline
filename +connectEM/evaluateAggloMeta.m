function y = evaluateAggloMeta(skelpath, graph, segmentMeta, agglos, p, nameOfRun, mainFolder, limitaggloNum, limitaggloSize, agglos_reverse, maxTube)

segmentMeta.point = segmentMeta.point';
doc_folder = [mainFolder, nameOfRun, '/'];
mkdir(doc_folder);
for file_idx = 1 : length(skelpath)
    % Load skeleton and shift to right coordinate system if traced in
    % another dataset
    idx = 1;
    % skel = skeleton(skelpath{file_idx});
    % assert(length(skel.nodes) == 1);
    % if isempty(strfind(skel.parameters.experiment.name, 'ROI2017'))
    %     skel.nodes{idx}(:, 1 : 3) = bsxfun(@minus, skel.nodes{idx}(:, 1 : 3),  [1195, 1515, 115]-(129 - [25 25 10]));
    % end
    % skel.nodes{idx}(skel.nodes{idx} <= 0) = 1;
    % skel.nodesNumDataAll{idx}(:, 3 : 5)= skel.nodes{idx}(:, 1 : 3);
    % % Make sure skeleton is one CC and has no circles
    % assert(graphconncomp(skel.createAdjacencyMatrix(idx)) == 1);
    % assert(~skel.hasCircle(idx));
    % % Look up segment IDs, remove all nodes iteratively that hit background
    % % and are of degree one (to remove parts outside segmented area)
    % ids = Seg.Global.getSegIds(p, skel.nodes{idx}(:, 1 : 3));
    % skel = connectEM.evaluateAggloCleanSkel(skel, idx, ids == 0);
    % % Could modify ids directly instead of looking up again
    % ids = Seg.Global.getSegIds(p, skel.nodes{idx}(:, 1 : 3));
    % %Main function calculating the metrics
    s = strsplit(skelpath{file_idx}, '/');
    load(['/gaba/scratch/kboerg/' s{end} '.mat'], 'skel', 'ids');
    [recall, splits, mergers, validnodes, foundAgglomerates, connM] = connectEM.evaluateAgglo(agglos, segmentMeta, skel, ids, graph.neighbours, limitaggloNum, limitaggloSize, agglos_reverse, maxTube);
    % Collect results in structure
    y.recall_col{file_idx} = recall;
    y.splits_col(file_idx) = splits;
    y.mergers_col(file_idx) = mergers;
    y.length_col(file_idx) = full(sum(sum(skel.createWeightedAdjacencyMatrix(idx)))) / 2;
    good_edges = skel.edges{idx}(all(ismember(skel.edges{idx}, validnodes), 2), :);
    y.covered_col(file_idx) = sum(sqrt(sum(bsxfun(@times, skel.nodes{idx}(good_edges(:, 1), 1 : 3) - skel.nodes{idx}(good_edges(:, 2), 1 : 3), [11.24, 11.24, 28]).^2, 2)));
    y.foundAgglomerates_col{file_idx} = foundAgglomerates;
    y.connM{file_idx} = connM;
    % Write output results
    %connectEM.skeletonFromAgglo(graph.edges, segmentMeta, agglos(foundAgglomerates), [num2str(file_idx) '_'], doc_folder)
    %skel.write([doc_folder, 'skel_' num2str(file_idx)  '.nml']);
end

end
