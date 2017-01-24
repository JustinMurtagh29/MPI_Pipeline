function [ interfaces, intIdx ] = calculateInterfaces( seg, edges, borders, areaT, voxelSize, rinclude )
%CALCULATEINTERFACES Calculate interfaces for synapse detection densely for
% a cube.
% INPUT seg: 3d int
%           Segmentation matrix.
%       edges: [Nx2] int
%           Edge list where edge row specifies the ids of two neighboring
%           segments.
%           (see output of SynEm.Svg.findEdgesAndBorders)
%       borders: [Nx3] table
%           Table with border information containing information about the
%           corresponding row in edges.
%           (see output of SynEm.Svg.findEdgesAndBorders)
%       areaT: int
%           Lower area threshold on border size. Only borders with
%           size > areaT are considered.
%       voxelSize: [1x3] double
%           Voxel size per pixel in nm.
%       rinclude: [1xN] double
%           Size of the different subsegments in nm.
% OUTPUT interfaces: [Nx1] struct
%           Voxel data of each interface. Fields are
%           surface: [Nx1] cell array where each cell contains the linear
%               indices of voxels w.r.t. seg of an interface surface as an
%               [Nx1] int array.
%           subseg: [1xN] cell array where N = length(rinclude). Each cell
%               contains a [Mx2] cell array where
%               M = length(interfaceSurface) which contains the linear
%               indices of the first and second subsegment w.r.t. seg.
%        intIdx: [Nx1] logical
%           Logical indices of the borders that were considered for
%           interface calculation.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

fprintf(['[%s] SynEM.Svg.calculateInterfaces - Starting interface ', ...
    'calculation.\n'], datestr(now));

%get borders above area threshold
area = [borders(:).Area];
intIdx = area > areaT;
interfaces.surface = {borders(intIdx).PixelIdxList}';
if isrow(interfaces.surface{1}) %always use column vectors
    interfaces.surface = cellfun(@(x)x',interfaces.surface, ...
        'UniformOutput',false);
end
neighborIDs =  edges(intIdx,:);

%calculate subsegments
interfaces.subseg = SynEM.Svg.calculateSubsegments( ...
    interfaces.surface, neighborIDs, seg, voxelSize, rinclude );
interfaces.rinclude = rinclude;

fprintf(['[%s] SynEM.Svg.calculateInterfaces - Finished interface ', ...
    'calculation.\n'], datestr(now));
end
