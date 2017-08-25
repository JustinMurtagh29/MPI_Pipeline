% load needed data
load(fullfile(outputFolder,'superagglos.mat'), 'superagglos');
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat','p');
[graph, segmentMeta] = connectEM.loadAllSegmentationData(p);
borderMeta = load([p.saveFolder 'globalBorder.mat'], 'borderCoM');

% create bbox
options.border = [2000; -2000];
borderNm = repmat(options.border, 1, 3);
borderVoxel = round(bsxfun(@times, 1./p.raw.voxelSize, borderNm));
bboxSmall = p.bbox + borderVoxel';

% create fake center of dataset field in borderCoM list for correspondences
borderMeta.borderCoM(end+1,:)=round(mean(bboxSmall,2)');

edgesCol = [];

for idx = 1 : length(superagglos)
    if mod(idx,100) == 0
        idx
    end
    nodesHere = superagglos(idx).nodes;
    segIdsHere = nodesHere(:, 4);
    % find edge idxs of used hot edges
    [~, usededges] = ismember(graph.edges, sort(segIdsHere(superagglos(idx).edges),2),'rows');
    % find whether edge is outside of bbox
    borderidxs = graph.borderIdx(usededges>0);
    borderidxs(isnan(borderidxs))=length(borderMeta.borderCoM); % correspondences
    borderPositions = double(borderMeta.borderCoM(borderidxs,:));
    outsideBbox = ~(all(bsxfun(@gt, borderPositions, bboxSmall(:, 1)'), 2) & ...
        all(bsxfun(@lt, borderPositions, bboxSmall(:, 2)'), 2));
    % remove all insulting edges
    todelete = graph.prob(usededges>0)<0.98&outsideBbox;
    usededges(usededges==0) = [];
    edgesHere = superagglos(idx).edges ...
        + sum(arrayfun(@(x)size(x.nodes,1),superagglos(1:idx-1))); % offset so that all nodes and edges go into one struct
    edgesHere(usededges(todelete), :) = [];
    edgesCol = [edgesCol; edgesHere];
end
nodesCol = cat(1, superagglos(:).nodes);
C = Graph.findConnectedComponents(edgesCol);
superagglosBorderSplit = cellfun(@(x)struct('nodes',nodesCol(x, :), 'edges', edgesCol(all(ismember(edgesCol,x),2),:)), C);

% add lonely nodes
superagglosBorderSplit = [superagglosBorderSplit; struct('nodes', arrayfun(@(x){nodesCol(x,:)},setdiff(1:size(nodesCol,1), cell2mat(C))),1)];
