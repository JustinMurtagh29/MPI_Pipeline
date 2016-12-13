function [ edges, borders, segments, borderAdjacency ] = ...
    findEdgesAndBorders( seg, excludeVoxels )
%FINDEDGESANDBORDER Calculate edges of the supervoxel graph and the
%corresponding border pixels based on 26 connectivity.
% INPUT seg: 3d int
%           3d array of integer containing a label/segmentation matrix.
%       excludeVoxels: (optional) 3d logical
%           Array of same size as seg which specified voxels which are
%           excluded as boundary voxels. Voxels that should be excluded
%           should be set to true.
%           (E.g. if only a partial segmentation is available the region
%           outside the segmentation can be masked).
%           (Default: all possible voxels are considered).
% OUTPUT edges: [Nx2] int
%           Edge list where edge row specifies the ids of two neighboring
%           segments.
%        borders: [Nx3] struct
%           Voxel data of the corresponding row in edges. Fields are
%           'PixelIdxList': The linear indices of all boundary voxels of
%               the respective edge. A voxel is considered on the boundary
%               if it does not belong to a segment (i.e. has value 0 in
%               seg) and has exactly two neighboring segment ids in its 26
%               neighborhod.
%           'Area': Double specifying  area of the corresponding border in
%               number of voxels.
%           'Centroid': [1x3] vector of double containing the centroid of
%               the border in coordinates of seg.
%        segments: [Nx1] struct
%           Voxel data of the segments in seg. Fields are
%           'Id': Segmentation id of the segment.
%           'PixelIdxList': Linear indices of all voxles of the respective
%               segment.
%        borderAdjacency: [Nx2] int
%           Border adjacency edge list. Each row specifies the linear
%           indices of two adjacent borders.
% Author: Manuel Berning <manuel.berning@brain.mpg.de>
% Modified by: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

% Pad array with a 1 voxel surround with a new unique value
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
nSegId = seg(nInd);
nSegId = sort(nSegId,1);

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
edges = cat(2,reshape(nSegId(lSegId(:,toKeep)),2,[]),cell2mat(pairs))';
indExtend = repelem(ind(toExtend), numComb);
ind = cat(1,ind(toKeep),indExtend);
[uid,~,c] = unique(indExtend);
%these indices belong to several edges
overlapInd = [zeros(sum(toKeep),1); c];

%get the edge for each ind (use sorting and grouping to avoid find)
[edges,se] = sortrows(edges);
ind = ind(se);
overlapInd = overlapInd(se);
diffEdge = find([true; any(diff(edges),2)]);
if diffEdge(end) < size(ind,1)
    diffEdge(end + 1) = size(ind,1) + 1;
end

%make edges contiguous
groupIdx = zeros(length(ind),1); %resulting edge index for each ind
%linear indices of overlapping borders
borderAdjacency = cell(length(uid),1);
currGroupId = 1;
edgesNew = zeros(3*size(diffEdge,1),2,'like',edges);
borderPixelMat = false(size(seg)); %avoid slow ismember for adjacency_list
fprintf('Border calculation:  0%%');
for i = 1:length(diffEdge) - 1
    
    %get all voxels of unique edge
    borderIdx = ind(diffEdge(i):(diffEdge(i+1) - 1));
    borderPixelMat(borderIdx) = true;
    borderOverlapIdx = overlapInd(diffEdge(i):(diffEdge(i+1) - 1));
    
    %construct adjacency between voxels (using 26 neighborhood)
    adjacency_list = bsxfun(@plus, borderIdx, vec)';
    adjacency_list = num2cell(adjacency_list, 1);
    adjacency_list = cellfun(@(x)x(borderPixelMat(x))', adjacency_list, ...
        'UniformOutput', false);
    borderPixelMat(borderIdx) = false;
    
    %split into connected components
    bordersForThisEdge = bfs(borderIdx,adjacency_list); %see below
    for j = 1:length(bordersForThisEdge)
        
        %add indices to border
        currGroup = ismember(borderIdx, bordersForThisEdge{j});
        currGroupIdx = diffEdge(i):(diffEdge(i+1) - 1);
        currGroupIdx = currGroupIdx(currGroup);
        groupIdx(currGroupIdx) = currGroupId;
        edgesNew(currGroupId,1:2) = edges(diffEdge(i),:);
        
        %get border adjacency edges
        if nargout == 4
            tmp = setdiff(borderOverlapIdx(currGroup),0);
            for k = 1:length(tmp)
                borderAdjacency{tmp(k)} = [borderAdjacency{tmp(k)}; ...
                                           currGroupId];
            end
        end
        
        currGroupId = currGroupId + 1;
    end
    if floor(i/(length(diffEdge) - 1)*100) < 10
        fprintf('\b\b%d%%',floor(i/(length(diffEdge) - 1)*100));
    else
        fprintf('\b\b\b%d%%',floor(i/(length(diffEdge) - 1)*100));
    end
end
fprintf('\n');

%final edges list
edgesNew(currGroupId:end,:) = [];
edges = edgesNew;

%final border adjacency list
if nargout == 4
    borderAdjacency(cellfun(@length,borderAdjacency) < 2) = [];
    borderAdjacency = cellfun(@(x)combnk(x,2),borderAdjacency, 'uni', 0);
    borderAdjacency = cell2mat(borderAdjacency);
end

%save corresponding wall voxels to borders
[groupIdx,I] = sort(groupIdx);
ind = ind(I);
groupIdx = [1; diff(groupIdx); 1];
groupIdx = find(groupIdx);
groupIdx = diff(groupIdx);
PixelIdxList = mat2cell(ind,groupIdx);

%calculate area and centroid (and remove padding in calcCentroid)
Area = cellfun(@length,PixelIdxList);
sizSeg = size(seg);
[PixelIdxList, Centroid] = cellfun(@(x)calcCentroid(x,sizSeg), ...
    PixelIdxList,'UniformOutput',false);

%save to output
borders = table(PixelIdxList,Area,Centroid);
borders = table2struct(borders); %compatibility (would prefer table....)

%calculate segments if required
if nargout > 2
    seg = seg(2:end-1,2:end-1,2:end-1); %remove the padding
    segments = regionprops(seg, seg, 'PixelIdxList', 'MinIntensity');
    segments(arrayfun(@(x)isempty(x.PixelIdxList),segments)) = [];
    [segments.Id] = segments.MinIntensity;
    segments = rmfield(segments,'MinIntensity');
end
end

function [tInd, c] = calcCentroid(ind, siz)
%Calculate the centroid for a list of linear indices.
% INPUT ind: [Nx1] array of integer containing linear indices.
%       siz: [1x3] array of integer containing the size of the cube which
%           the linear indices refer to.
% OUTPUT tInd: [Nx1] array of integer as same length as ind containing the
%           linear indices in the coordinates of seg (without the padding).
%        c: [1x3] array of double containing the center of mass/mean
%           coordinates as x-, y- and z-coordinate.

[x,y,z] = ind2sub(siz,ind);
c = mean([x,y,z] - 1); %remove padding
tInd = sub2ind(siz - 2,x - 1,y - 1,z - 1); %remove padding
end

function components = bfs( idx, adjList )
%BFS Breadth-first search for boundary connected component calculation.
% INPUT idx: [Nx1] array of integer containing a list of linear indices
%           for which the connected components are calculated.
%       adjacency_list: [Nx1] cell array of integer containing the adjacend
%           linear indices for each index in idx as a row vector, i.e. the
%           i-th cell contains all linear indices which are adjacent to
%           idx(i) as an [1xN] array.
% OUTPUT components: [Nx1] cell array. Each cell contains the linear
%           indices of one connected component within idx.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

components = cell(0,1);
k = 1;
nonvisited = true(length(idx), 1);
idxLookup = sparse(double(idx), 1, 1:length(idx));
currComp = false(length(idx), 1);
queuedNodes = false(length(idx), 1);

while any(nonvisited)
    currComp = currComp & false;
    
    %get first node of components
    nextIdx = find(nonvisited, 1, 'first');
    currComp(nextIdx) = true;
    queuedNodes(nextIdx) = true;
    
    while any(queuedNodes)
        queuedNodes(idxLookup(horzcat(adjList{queuedNodes}))) = true;
        queuedNodes = queuedNodes & ~currComp;
        currComp = currComp | queuedNodes;
    end
    components{k, 1} = idx(currComp)';
    nonvisited(currComp) = false;
    k = k + 1;
end

end


