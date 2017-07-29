function [ wallIndices,uneighboursID] = getContactIndices( segments,walls )
% GETCONTACTIDX finds all contact voxels that have exactly two neighboring
% segments and returns them together with the IDs of those segments.
%
% OUTPUT: wallIndices: Array containing the wall voxels.
%         uneighboursID: 2xn array where each row contains the segment IDs
%                        for the neighbors of the corresponding wall voxel
%                        in wallIndices.

cubesize = size(segments);
nD = ndims(segments);

fprintf('%s INFO calculating contacts and neighboring cell IDs\n',datestr(now,'yyyy-mm-dd HH:MM:SS,FFF'));
 
% clear all to close to border pixel
walls(1,:,:) = false;
walls(end,:,:) = false;
walls(:,1,:) = false;
walls(:,end,:) = false;
walls(:,:,1) = false;
walls(:,:,end) = false;
wallIndices = find(walls);
 
% generate indexshifters for 26 neighborhood
[x,y,z] = meshgrid([-1 0 1], [-1 0 1]*cubesize(1), [-1 0 1]*prod(cubesize(1:end-1))); %linear indices for neighboring voxels in all three dimensions
n26 = x+y+z;
n26 = n26([1:13, 15:end]);
 
neighbourIndices = bsxfun(@plus,wallIndices,n26); % calculate indices of neighbouring voxels for all wall voxels %rows of wallidxi are the linear indices of the neighboring voxels
neighbourIDs = sort(segments(neighbourIndices), 2); % sort neighboring  segments IDs for all wall wall voxels %rows of neighboursID contain the ids of the segments of the voxels sorted in ascending order
uneighbours = [true(size(neighbourIDs,1),1), 0<diff(neighbourIDs,1,2)]; % and find the unique neighboring segments IDs for all wall voxels %rows of uneighbours contain 0 or 1 for whether the id increases and 1 at the beginning
leftout = sum(uneighbours,2) ~= 3; % find junctions (<10% ) %take only rows which have exactly two jumps (plus the additional 1 at the beginning of each row)
wallIndices = wallIndices(~leftout);% forget junctions (<10% )
uneighbours = uneighbours(~leftout,:);
neighbourIDs = neighbourIDs(~leftout,:).';
uneighboursID = reshape(neighbourIDs(uneighbours.'), nD, []).';
 
uneighboursID = uneighboursID(:,2:3);% forget wall id
end