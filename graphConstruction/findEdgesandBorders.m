function [seg, edges, borders, edgesToBorder] = findEdgesandBorders(seg)
% Computation of edges and borders optimized for 512x512x256

% Pad array with a 1 voxel surround with a new unique value
globalBorderId = max(seg(:))+1;
seg = padarray(seg,[1 1 1],globalBorderId);

% Size of padded segmentation
[M,N,P] = size(seg);

% Construct 26-connectivity linear indices shift for padded segmentation
vec = int32([(-M*N+[-M-1 -M -M+1 -1 0 1 M-1 M M+1]) [-M-1 -M -M+1 -1 1 M-1 M M+1] (M*N+[-M-1 -M -M+1 -1 0 1 M-1 M M+1])]);

% Find linear inidices of all wall voxel
ind = int32(find(seg==0));

% Find segmentation ID of all neighbours of all wall voxel (according to 26
% connectivity)
nInd = bsxfun(@plus, ind, vec);
nSegId = seg(nInd);

% Find edges
edges = findEdges(nSegId);

% Remove edges to padded value (see first comment)
% This uses fact that globalBorderId is largest as is made sure in the beginning 
edgesToBorder = edges(edges(:,2) == globalBorderId,:);
edges(edges(:,2) == globalBorderId,:) = [];

%Find leaves
leaves = findLeaves(edges,edgesToBorder(:,1));

% Find borders
[edges,borders] = findBorders(ind, nInd, nSegId, edges, M, N, P);

%Undo padding of seg
seg = seg(2:end-1,2:end-1,2:end-1);

%Delete leaves from seg, edges and borders
[seg,edges,borders] = correctLeaves(seg,leaves,edges,borders);

end
