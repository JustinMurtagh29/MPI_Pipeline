function findEdgesandBorders(segFile, edgeFile, borderFile, segmentFile, tileBorder)
% Computation of edges and borders optimized for 512x512x256

% Load segmentation from file
[seg, segSmall] = loadSegData(segFile, tileBorder);

% Pad array with a 1 voxel surround with a new unique value
globalBorderId = max(seg(:))+1;
segSmall = padarray(segSmall,[1 1 1],globalBorderId);

% Size of padded segmentation
[M,N,P] = size(segSmall);

% Construct 26-connectivity linear indices shift for padded segmentation
vec = int32([(-M*N+[-M-1 -M -M+1 -1 0 1 M-1 M M+1]) [-M-1 -M -M+1 -1 1 M-1 M M+1] (M*N+[-M-1 -M -M+1 -1 0 1 M-1 M M+1])]);

% Find linear inidices of all wall voxel
ind = int32(find(segSmall==0));

% Find segmentation ID of all neighbours of all wall voxel (according to 26
% connectivity)
nInd = bsxfun(@plus, ind, vec);
nSegId = segSmall(nInd);

% Find edges
edges = findEdges(nSegId);

% Remove edges to padded value (see first comment)
% This uses fact that globalBorderId is largest as is made sure in the beginning 
edgesToBorder = edges(edges(:,2) == globalBorderId,:);
edges(edges(:,2) == globalBorderId,:) = [];

% Find leaves
leaves = findLeaves(edges,edgesToBorder(:,1));

% Find borders
[edges,borders] = findBorders(ind, nInd, nSegId, edges, M, N, P);

%Undo padding of seg
segSmall = segSmall(2:end-1,2:end-1,2:end-1);

%Delete leaves from seg, edges and borders
[segSmall,edges,borders] = correctLeaves(segSmall,leaves,edges,borders);

% Save to files: currently segmentation overwrites old one (now leaves are merged)
seg(1-tileBorder(1,1):end-tileBorder(1,2),...
	1-tileBorder(2,1):end-tileBorder(2,2),...
	1-tileBorder(3,1):end-tileBorder(3,2)) = segSmall;
save(segFile, 'seg');
save(edgeFile, 'edges', 'edgesToBorder');
save(borderFile, 'borders');

% get segment PixelIdxLists for glia predition
ids = unique(segSmall);
props = regionprops(segSmall,'PixelIdxList');
ids(ids == 0) = [];
for i = 1:length(ids)
    segments(i).PixelIdxList = props(ids(i)).PixelIdxList;
    segments(i).id = ids(i);
end
save(segmentFile,'segments','-v7.3');

end

