function [edges, ind] = codeBenedikt(seg)
borderId = max(seg(:)) + 1;
seg = padarray(seg,[1, 1, 1],borderId);
[M,N,~] = size(seg);

% Construct 26-connectivity linear indices shift for padded segmentation
vec = int32([(-M*N+[-M-1 -M -M+1 -1 0 1 M-1 M M+1]) ...
    [-M-1 -M -M+1 -1 1 M-1 M M+1] (M*N+[-M-1 -M -M+1 -1 0 1 M-1 M M+1])]);

% Find linear inidices of all wall voxel
if exist('excludeVoxels','var') && ~isempty(excludeVoxels)
    excludeVoxels = padarray(excludeVoxels,[1, 1, 1],0);
    ind = int32(find(seg == 0 & ~excludeVoxels));
else
    ind = int32(find(seg==0));
end

% Find segmentation ID of all neighbours of all wall voxel (according to 26
% connectivity)
nInd = bsxfun(@plus, ind', vec');
seg(padarray(false(size(seg) - 2),[1, 1, 1],1)) = 0; % set boundary to zero
nSegId = seg(nInd);  % get segment IDs of neighbors
nSegId = sort(nSegId,1); % sort segment IDs along 1st dim

%get neighbors for wall voxels
lSegId = [false(1,size(nSegId,2)); diff(nSegId,1,1)>0];
numNeighbors = sum(lSegId,1);
toKeep = (numNeighbors == 2)';
toExtend = (numNeighbors > 2)';

%get combinations for voxels with between more than 2 ids
nSegIdExtend = num2cell(nSegId(:,toExtend),1);
nSegId = nSegId(:,toKeep);
nSegIdExtend = cellfun(@(x,ind)x(ind),nSegIdExtend, ...
    num2cell(lSegId(:,toExtend),1), 'UniformOutput', false);
pairs = cellfun(@(x)combnk(x,2)',nSegIdExtend,'UniformOutput',false);
numComb = cellfun(@(x)size(x,2),pairs)';

%save edge for each voxel and voxel linear indices (ind)
edges = cat(2,reshape(nSegId(lSegId(:,toKeep)),2,[]),cell2mat(pairs))';  %reshape nSegIds with only two neighbors and add pair list
indExtend = repelem(ind(toExtend), numComb,1);
ind = cat(1,ind(toKeep), indExtend);
[uid,~,c] = unique(indExtend);
%these indices belong to several edges
overlapInd = [zeros(sum(toKeep),1); c];  % 0 where ind belongs to one edge and x where ind is uid(x)
