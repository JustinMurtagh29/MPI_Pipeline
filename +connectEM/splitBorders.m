% load needed data
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat','p');
temp = load(fullfile(p.saveFolder,'aggloState/axons_03.mat'));
superagglos = temp.axons(temp.indBigAxons);
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
    probs = cell2mat(graph.neighProb(segIdsHere));
    % find whether edge is outside of bbox
    borderidxs = cell2mat(graph.neighBorderIdx(segIdsHere));
    borderidxs(isnan(borderidxs))=length(borderMeta.borderCoM); % correspondences
    borderPositions = double(borderMeta.borderCoM(borderidxs,:));
    outsideBbox = ~(all(bsxfun(@gt, borderPositions, bboxSmall(:, 1)'), 2) & ...
        all(bsxfun(@lt, borderPositions, bboxSmall(:, 2)'), 2));
    % remove all insulting edges
    actualEdges = [repelem(segIdsHere, cellfun('length',graph.neighbours(segIdsHere))), cell2mat(graph.neighbours(segIdsHere))];
    actualEdges(probs<0.98&outsideBbox,:) = [];
    edgesHere = actualEdges ...
        + sum(arrayfun(@(x)size(x.nodes,1),superagglos(1:idx-1))); % offset so that all nodes and edges go into one struct
    edgesCol = [edgesCol; edgesHere];
end
nodesCol = cat(1, superagglos(:).nodes);
C = Graph.findConnectedComponents(edgesCol);
superagglosBorderSplit = cellfun(@(x)struct('nodes',nodesCol(x, :), 'edges', edgesCol(all(ismember(edgesCol,x),2),:)), C);

% add lonely nodes
superagglosBorderSplit = [superagglosBorderSplit; struct('nodes', arrayfun(@(x){nodesCol(x,:)},setdiff(1:size(nodesCol,1), cell2mat(C))),1)];

calculateLength = @(x)max(pdist(bsxfun(@times, double(borderMeta.borderCoM(x, :)), p.raw.voxelSize)));
filternan = @(x)x(~isnan(x));
for idx = 1 : length(superagglosBorderSplit)
    idx
    axonLength(idx) = max([-1, calculateLength(filternan(cell2mat(graph.neighBorderIdx(superagglosBorderSplit(idx).nodes(:,4)))))]);
end
