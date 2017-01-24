function [ data, metadata ] = convertForSynapseAnnotationGUI( raw, seg, ...
    edges, interfaces, intRefSiz )
%CONVERTFORSYNAPSEANNOTATIONGUI Convert interfaces for synapse annotation
%GUI.
% INPUT raw: 3d int
%           Raw data images. Should have a border of [100, 100, 40] on each
%           side compared to the cube on which the interfaces were
%           calculated. Otherwise raw will be padded with zeros.
%       seg: 3d int
%           Volume segmentation for raw data. Should have a border of
%           [100, 100, 40] on each side compared to the cube on which the
%           interfaces were calculated. Otherwise seg will be padded with
%           zeros.
%       edges: [Nx2] int
%           see SynEM.Svg.findEdgesAndBorders
%       interfaces: struct
%           see SynEM.Svg.calculateInterfaces
%       intRefSiz: (Optional) [1x3] int
%           Size of the segmentation matrix for which the interfaces were
%           calculated.
%           (Default: [512 512 256])
% OUTPUT data: struct
%           Data struct required by SynapseAnnotationGUI.
%        metadata: struct
%           Metadata struct containing additional information about data.
%           (Note that this function leaves most fields empty which are not
%           required for the GUI but contain additional information about
%           the data struct.)
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('intRefSiz', 'var') || isempty(intRefSiz)
    intRefSiz = [512 512 256];
end
if any((size(raw) - intRefSiz) ~= [200, 200, 80])
    raw = padarray(raw,[100, 100, 40]);
end
if any((size(seg) - intRefSiz) ~= [200, 200, 80])
    seg = padarray(seg,[100, 100, 40]);
end

data.raw = raw;
data.segments = seg;
data.interfaceSurfaceList = interfaces.surface;
data.subsegmentsList = interfaces.subseg;

%convert indices stored in interfaces to larger cube
for i = 1:length(data.interfaceSurfaceList)
    [x, y, z] = ind2sub(intRefSiz, data.interfaceSurfaceList{i});
    x = x + 100;
    y = y + 100;
    z = z + 40;
    data.interfaceSurfaceList{i} = uint32(sub2ind(size(raw), x, y, z));
    for j = 1:2
        [x, y, z] = ind2sub(intRefSiz, data.subsegmentsList{1}{i,j});
        x = x + 100;
        y = y + 100;
        z = z + 40;
        data.subsegmentsList{1}{i,j} = uint32(sub2ind(size(raw), x, y, z));
    end
end
data.neighborIDs = edges;

%initialize metadata fields
metadata.experiment = [];
metadata.segmentation = [];
metadata.boundingBox = [];
metadata.innerBoundingBox = [];
metadata.rinclude = interfaces.rinclude;

end

