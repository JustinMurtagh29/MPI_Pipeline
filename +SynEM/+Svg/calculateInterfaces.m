function [ interfaces, intIdx, edges, borders ] = calculateInterfaces( ...
    seg, edges, borders, areaT, voxelSize, rinclude )
%CALCULATEINTERFACES Calculate interfaces for synapse detection densely for
% a cube.
% INPUT seg: 3d int
%           Segmentation matrix.
%       edges: (Optional) [Nx2] int
%           Edge list where edge row specifies the ids of two neighboring
%           segments.
%           (Default: will be calculated from seg using
%                     SynEm.Svg.findEdgesAndBorders)
%       borders: (Optional) [Nx3] table
%           Table with border information containing information about the
%           corresponding row in edges.
%           (Default: will be calculated from seg using
%                     SynEm.Svg.findEdgesAndBorders)
%       areaT: int
%           Lower area/size threshold on border size. Only borders with
%           size > areaT are considered.
%       voxelSize: [1x3] double
%           Voxel size per pixel in nm.
%       rinclude: [1xN] double
%           Size of the different subsegments in nm.
% OUTPUT interfaces: struct
%           Voxel data of each interface. Fields are
%           surface: [Nx1] cell array where each cell contains the linear
%               indices of voxels w.r.t. seg of an interface surface as an
%               [Nx1] int array.
%           subseg: [1xN] cell array where N = length(rinclude). Each cell
%               contains a [Mx2] cell array where
%               M = length(interfaceSurface) which contains the linear
%               indices of the first and second subsegment w.r.t. seg.
%        intIdx: [Nx1] logical
%           Logical indices of the edges and borders that were considered
%           for interface calculation after applying the size threshold.
%        edges: see input
%        bordes: see input
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

fprintf(['[%s] SynEM.Svg.calculateInterfaces - Starting interface ', ...
    'calculation.\n'], datestr(now));

if isempty(edges) || isempty(borders)
    fprintf(['[%s] SynEM.Svg.calculateInterfaces - Calculating ', ...
        'edges and borders.\n'], datestr(now));
    [edges, borders] = SynEM.Svg.findEdgesAndBorders(seg);
end

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
