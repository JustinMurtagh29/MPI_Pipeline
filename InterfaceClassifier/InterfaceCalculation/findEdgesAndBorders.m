function [ edges, borders ] = findEdgesAndBorders( seg, excludeVoxels, mergeMode )
%FINDEDGESANDBORDER Calculate edges of the supervoxel graph and the
%corresponding border pixels.
% INPUT seg: 3d array of integer containing a label/segmentation matrix.
%       excludeVoxels: (Optional) 3d array of logical of same size as seg
%           which specified voxels which are excluded as boundary voxels.
%           Voxels that should be excluded should be set to true.
%           (E.g. if only a partial segmentation is available the region
%           outside the segmentation can be maske).
%           (Default: all possible voxels are considered).
%       mergeMode: (Optional) Logical specifying whether wall voxels with
%           exactly one segment ID in their 26 neighborhood should be added
%           to this ID.
%           (Default: false)
% OUTPUT edges: [Nx2] array of integer where each row contains an
%           supervoxel graph edge.
%        borders: [Nx3] table containing the rows.
%           'PixelIdxList': The linear indices of all boundary voxels of
%               the respective edge. A voxel is considered on the boundary
%               if it does not belong to a segment (i.e. has value 0 in
%               seg) and has exactly two neighboring segment ids in its 26
%               neighborhod.
%           'Area': Double specifying  area of the corresponding border in
%               number of voxels.
%           'Centroid': [1x3] vector of double containing the centroid of
%               the border in coordinates of seg.
% Author: Manuel Berning <manuel.berning@brain.mpg.de>
% Modified by: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('mergeMode','var') || isempty(mergeMode)
    mergeMode = false;
end

% Pad array with a 1 voxel surround with a new unique value
borderId = max(seg(:)) + 1;
seg = padarray(seg,[1, 1, 1],borderId);
[M,N,~] = size(seg);

% Construct 26-connectivity linear indices shift for padded segmentation
vec = int32([(-M*N+[-M-1 -M -M+1 -1 0 1 M-1 M M+1]) [-M-1 -M -M+1 -1 1 M-1 M M+1] (M*N+[-M-1 -M -M+1 -1 0 1 M-1 M M+1])]);

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

if mergeMode
    %set voxels with exactly one neighbor to this neighbor ID
    nSegId = seg(nInd);
    nSegId = sort(nSegId,1);
    lSegId = [false(1,size(nSegId,2)); diff(nSegId,1,1)>0];
    toMerge = sum(lSegId,1) == 1;
    
    toMergeId = nSegId(:,toMerge);
    toMergeId = toMergeId(lSegId(:,toMerge));
    toMerge(toMerge) = toMergeId ~= borderId;
    
    toMergeId = nSegId(:,toMerge);
    toMergeId = toMergeId(lSegId(:,toMerge));
    toMergeInd = ind(toMerge);
    
    nSegId(:,toMerge) = [];
    nSegId(nSegId == borderId) = 0;
    ind(toMerge) = [];

    seg(toMergeInd) = toMergeId;
    nSegId(:,toMerge) = [];
else
    seg(padarray(false(size(seg) - 2),[1, 1, 1],1)) = 0; %set boundary to zero again
    nSegId = seg(nInd);
    nSegId = sort(nSegId,1);
end

%get wall voxels with exactly two neighboring segment ids
lSegId = [false(1,size(nSegId,2)); diff(nSegId,1,1)>0];
keep = sum(lSegId,1) == 2;

%delete other wall voxel indices
ind = ind(keep);
nSegId = nSegId(:,keep);
lSegId = lSegId(:,keep);

%get the edge for each ind
edges = reshape(nSegId(lSegId), 2, size(lSegId,2))';
[edges,se] = sortrows(edges);
ind = ind(se);
diffEdge = find([true; any(diff(edges),2)]);
if diffEdge(end) < size(ind,1)
    diffEdge(end + 1) = size(ind,1) + 1;
end

%make edges contiguous
groupIdx = zeros(length(ind),1); %resulting edge index for each ind
currGroupId = 1;
edgesNew = zeros(3*size(diffEdge,1),2,'like',edges);
for i = 1:length(diffEdge) - 1
    %get all voxels of unique edge
    borderIdx = ind(diffEdge(i):(diffEdge(i+1) - 1));
    %construct adjacency between voxels (using 26 neighborhood)
    adjacency_list = bsxfun(@plus, borderIdx, vec)';
    adjacency_list = mat2cell(adjacency_list, 26, ones(size(adjacency_list,2),1));
    adjacency_list = cellfun(@(x)x(ismember(x,borderIdx))', adjacency_list, 'UniformOutput', false);
    %split into connected components
    bordersForThisEdge = bfs(borderIdx,adjacency_list);
    for j = 1:length(bordersForThisEdge)
        currGroup = ismember(borderIdx, bordersForThisEdge{j});
        currGroupIdx = diffEdge(i):(diffEdge(i+1) - 1);
        currGroupIdx = currGroupIdx(currGroup);
        groupIdx(currGroupIdx) = currGroupId;
        edgesNew(currGroupId,1:2) = edges(diffEdge(i),:);
        currGroupId = currGroupId + 1;
    end
end
edgesNew(currGroupId:end,:) = [];
edges = edgesNew;

%save corresponding wall voxels to borders
[groupIdx,I] = sort(groupIdx);
ind = ind(I);
groupIdx = [1; diff(groupIdx); 1];
groupIdx = find(groupIdx);
groupIdx = diff(groupIdx);
PixelIdxList = mat2cell(ind,groupIdx);

%calculate area and centroid
Area = cellfun(@length,PixelIdxList);
sizSeg = size(seg);
[PixelIdxList, Centroid] = cellfun(@(x)calcCentroid(x,sizSeg),PixelIdxList,'UniformOutput',false);

%save to output
borders = table(PixelIdxList,Area,Centroid);
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
c = mean([x,y,z]);
tInd = sub2ind(siz - 2,x - 1,y - 1,z - 1);
end

function components = bfs(idx, adjacency_list)
% Calculate the connected components using a breadth-first search.
% INPUT idx: [Nx1] array of integer containing a list of linear indices
%           for which the connected components are calculated.
%       adjacency_list: [Nx1] cell array of integer containing the adjacend
%           linear indices for each index in idx as a row vector, i.e. the
%           i-th cell contains all linear indices which are adjacent to
%           idx(i) as an [1xN] array.
% OUTPUT components: [Nx1] cell array. Each cell contains the linear
%           indices of one connected component within idx.

components ={};
k = 1;
nonvisited = 1:length(idx);

while ~isempty(nonvisited)
    components{k,1} = idx(nonvisited(1));
    l1 = size(components{k,1},2);
    components{k,1}=horzcat(components{k,1},adjacency_list{nonvisited(1)});
    l2 = size(components{k,1},2);
    nonvisited(1)=[];
    while l1~=l2
        v=[];
        for j=l1+1:l2
            index = find(idx==(components{k,1}(j)));
            v=[v,adjacency_list{index}];
            nonvisited(nonvisited==index)=[];
        end
        l1 = l2;
        v = intersect(v,idx(nonvisited))';
        components{k,1}=horzcat(components{k,1},v);
        l2 = size(components{k,1},2);
    end
    k = k+1;
end
end