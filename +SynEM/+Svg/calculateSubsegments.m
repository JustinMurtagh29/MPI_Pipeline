function subsegmentsList = calculateSubsegments( interfaceSurfaceList, neighborIDs, seg, voxelSize, rinclude  )
%CALCULATESUBSEGMENTS Calculate parts of segments adjacent to borders.
% INPUT interfaceSurfaceList: [Nx1] cell of [Mx1] int
%           Each cell contains all linear indices of all voxel from a
%           single border.
%       neighborIDs: [Nx2] int
%           Edge list for the borders in interfaceSurfaceList. Each row
%           contains the segment ids of the corresponding border in
%           interfaceSurfaceList.
%       seg: 3d int
%           Segmentation array.
%       voxelSize: [1x3] double
%           Voxel size per pixel in nm.
%       rinclude: [1xN] double
%           Size of the different subsegments in nm.
% OUTPUT subsegmentsList: [1xN] cell of [Mx2] cells of [Kx1] int
%           The outermost cell array corresponds to the different
%           subsegment sizes specified via rinclude.
%           For each value in rinclude the [Mx2] cell contains the linear
%           indices of the two subsegments.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~isrow(rinclude)
    rinclude = rinclude';
end
rinclude = sort(rinclude);
if ~isrow(voxelSize)
    voxelSize = voxelSize';
end

fprintf(['[%s] SynEM.Svg.calculateSubsegments - Starting subsegment ', ...
    'calculation.\n'], datestr(now));

maxDist = ceil(max(rinclude)./voxelSize);
segSize = size(seg);

%prepare dilation mask (applied sequentially over rinclude)
h = arrayfun(@(x)sphereAverage(x, voxelSize), ...
    [rinclude(1), diff(rinclude)], 'uni', 0);

%prepare output
subsegmentsList = arrayfun(@(x)cell(length(interfaceSurfaceList), 2), ...
    1:length(rinclude), 'UniformOutput', false);

%calculate subsegments for each interface
fprintf('Subsegment calculation:  0%%');
for i = 1:length(interfaceSurfaceList)

    %get contact wall subscript indices
    [x,y,z] = ind2sub(segSize,interfaceSurfaceList{i});
    wallSubIndices = [x,y,z];

    %discard voxels too far away from interface surface
    roughframe = [max(min(wallSubIndices,[],1) - maxDist,[1 1 1]); ...
                  min(max(wallSubIndices,[],1) + maxDist,segSize)];

    %translate everything to local cube of size roughframe
    localSurfaceMask = false(diff(roughframe) + 1);
    relativeInterfaceCoords = bsxfun(@minus,wallSubIndices, ...
        roughframe(1,:) - 1);
    relativeInterfaceIndices = sub2ind(size(localSurfaceMask), ...
        relativeInterfaceCoords(:,1),relativeInterfaceCoords(:,2), ...
        relativeInterfaceCoords(:,3));
    localSurfaceMask(relativeInterfaceIndices) = true;
    localSeg = seg(roughframe(1):roughframe(2),...
                   roughframe(3):roughframe(4), ...
                   roughframe(5):roughframe(6));

    %calculate subsegments on local cube and translate indices to large
    %cube
    for k = 1:length(rinclude)
        localSurfaceMask = imdilate(localSurfaceMask,h{k});
        for j = 1:2
            subsegRelativeIndices = find(localSurfaceMask & ...
                (localSeg == neighborIDs(i,j)));
            [x,y,z] = ind2sub(size(localSurfaceMask),subsegRelativeIndices);
            subsegRelativeCoords = [x,y,z];
            subsegCoords = bsxfun(@plus, subsegRelativeCoords, ...
                roughframe(1,:) - 1);
            subsegmentsList{k}{i,j} = uint32(sub2ind(segSize, ...
                subsegCoords(:,1),subsegCoords(:,2),subsegCoords(:,3)));
        end
    end
    
    %print progress
    if floor(i/length(interfaceSurfaceList)*100) < 10
        fprintf('\b\b%d%%',floor(i/length(interfaceSurfaceList)*100));
    else
        fprintf('\b\b\b%d%%',floor(i/length(interfaceSurfaceList)*100));
    end
end
fprintf('\n');

fprintf(['[%s] SynEM.Svg.calculateSubsegments - Finished subsegment ', ...
    'calculation.\n'], datestr(now));
end

function h = sphereAverage(radius, voxelSize)
%SPHEREAVERAGE Create a sphere shaped binary array.
% INPUT radius: Double specifying the sphere radius.
%       voxelSize: [3x1] vector of double specifying the voxel size in nm.
% OUTPUT h: 3d spherical mask.

siz = ceil(radius./voxelSize);
[x,y,z] = meshgrid(-siz(1):siz(1),-siz(2):siz(2),-siz(3):siz(3));
h = (voxelSize(1).*x).^2 + (voxelSize(2).*y).^2 + ...
    (voxelSize(3).*z).^2 <= (radius)^2;
end
