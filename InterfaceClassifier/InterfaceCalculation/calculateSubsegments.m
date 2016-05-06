function subsegmentsList = calculateSubsegments( interfaceSurfaceList, neighborIDs, seg, d, voxelSize )
%CALCULATESUBSEGMENTS Calculate parts of segments adjacent to borders.
% INPUT interfaceSurfaceList: [Nx1] cell array where each entry contains
%           the linear indices of a border/interface surface with respect
%           to seg.
%       neighborIDs: [Nx2] array of integer where each row contains the ids
%           of the neighboring segments of the corresponding interface in
%           interfaceSurfaceList.
%       seg: 3d array of integer containing the segmentation.
%       d: [1xN] array of integer containing the maximal distances of the
%           respective subsegment to the interface surface in nm.
%       voxelSize: [1x3] array containing the voxel size in nm.
% OUTPUT subsegmentsList: [1xN] cell array where N = length(d)
%           corresponding to the sizes of the subsegments. Each cell
%           contains another cell array of size
%           length(interfaceSurfaceList)x2 containing the linear indices of
%           the subsegments for the corresponding interface surface.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~isrow(d)
    d = d';
end
if ~isrow(voxelSize)
    voxelSize = voxelSize';
end

fprintf('[%s] Restricting segments to parts close to interface surface.\n',datestr(now));

maxDist = ceil(max(d)./voxelSize);
segSize = size(seg);
h = cell(length(d),1);

subsegmentsList = cell(1,length(d));
for k = 1:length(subsegmentsList)
    subsegmentsList{k} = cell(length(interfaceSurfaceList),2);
    h{k} = sphereAverage(d(k), voxelSize);
end

for iInterface = 1:length(interfaceSurfaceList)

    %get contact wall subscript indices
    [x,y,z] = ind2sub(segSize,interfaceSurfaceList{iInterface});
    wallSubIndices = [x,y,z];

    % every location out of this from cannot be within rinclude to any contact voxel
    roughframe = [max(min(wallSubIndices,[],1) - maxDist,[1 1 1]); min(max(wallSubIndices,[],1) + maxDist,segSize)];
    
    %translate everything to local cube of size roughframe
    iCube = false(diff(roughframe) + 1);
    relativeInterfaceCoords = wallSubIndices - repmat(roughframe(1,:),size(wallSubIndices,1),1) + 1;
    relativeInterfaceIndices = sub2ind(size(iCube),relativeInterfaceCoords(:,1),relativeInterfaceCoords(:,2),relativeInterfaceCoords(:,3));
    iCube(relativeInterfaceIndices) = true;
    localSeg = seg(roughframe(1):roughframe(2),roughframe(3):roughframe(4),roughframe(5):roughframe(6));
    
    %calculate subsegments on local cube and translate indices to large
    %cube
    for k = 1:length(d)
        dilatedContact = imdilate(iCube,h{k});
        for iCell = 1:2
            subsegRelativeIndices = find(dilatedContact & (localSeg == neighborIDs(iInterface,iCell)));
            [x,y,z] = ind2sub(size(dilatedContact),subsegRelativeIndices);
            subsegRelativeCoords = [x,y,z];
            subsegCoords = subsegRelativeCoords + repmat(roughframe(1,:),size(subsegRelativeCoords,1),1) - 1;
            subsegmentsList{k}{iInterface,iCell} = uint32(sub2ind(segSize,subsegCoords(:,1),subsegCoords(:,2),subsegCoords(:,3)));
        end
    end
end
end

function h = sphereAverage(radius, voxelSize)
%SPHEREAVERAGE Create a sphere shaped binary array.
% INPUT radius: Double specifying the sphere radius.
%       voxelSize: [3x1] vector of double specifying the voxel size in nm.
% OUTPUT h: 3d spherical mask.

siz = ceil(radius./voxelSize);
[x,y,z] = meshgrid(-siz(1):siz(1),-siz(2):siz(2),-siz(3):siz(3));
h = (voxelSize(1).*x).^2 + (voxelSize(2).*y).^2 + (voxelSize(3).*z).^2 <= (radius)^2;
end

