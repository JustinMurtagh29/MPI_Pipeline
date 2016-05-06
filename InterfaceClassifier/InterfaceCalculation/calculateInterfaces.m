function [ interfaces, intIdx ] = calculateInterfaces( seg, edges, borders, areaT, voxelSize, rinclude )
%CALCULATEINTERFACES Calculate interfaces for synapse detection.
% INPUT seg: 3d array of integer containing a segmentation.
%       edges: edge variable from the segmentation edge file.
%       borders: borders variable from the segmentation border file.
%       areaT: Integer specifying a lower area threshold on border size.
%           Only borders with size > areaT are considered.
%       voxelSize: [1x3] array of double containing the voxel size per
%           pixel in nm.
%       rinclude: [1xN] array of double containing the distance of the
% OUTPUT interfaces: Struct containing the fields
%           surface: [Nx1] cell array where each cell contains the linear
%               indices of voxels w.r.t. seg of an interface suface.
%           subseg: [1xN] cell array where N = length(rinclude). Each cell
%               contains a [Mx2] cell array where
%               M = length(interfaceSurface) which contains the linear
%               indices of the first and second subsegment w.r.t. seg.
%        intIdx: [Nx1] array of logical specifying the borders that were
%           considered for interface calculation.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%get borders above area threshold
area = [borders(:).Area];
intIdx = area > areaT;
interfaces.surface = {borders(intIdx).PixelIdxList}';
if isrow(interfaces.surface{1})
interfaces.surface = cellfun(@(x)x',interfaces.surface,'UniformOutput',false);
end
neighborIDs =  edges(intIdx,:);

%calculate subsegments
interfaces.subseg = calculateSubsegments(interfaces.surface, neighborIDs, seg, rinclude, voxelSize);
end
