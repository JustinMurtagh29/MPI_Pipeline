function [ data, metadata, interfaceLabels, undecidedList, voxelLabels ] = loadInterfaces( path, areaThreshold )
%LOADINTERFACES Load interface data file.
% INPUT path: File path.
%       areaThresold: Optional parameter specifying an area threshold to
%                     use. All interfaces greater than the threshold are
%                     kept. (Default: 150)

if ~exist('areaThreshold','var') || isempty(areaThreshold)
    areaThreshold = 150;
end

m = matfile(path);
data = m.data;
metadata = m.metadata;
interfaceLabels = m.interfaceLabels;
undecidedList = m.undecidedList;
try
    voxelLabels = m.voxelLabels;
catch
    voxelLabels = false(size(data.raw));
end

keepIndices = cellfun(@(x)length(x) > areaThreshold, data.interfaceSurfaceList);
data.interfaceSurfaceList = data.interfaceSurfaceList(keepIndices);
data.subsegmentsList = cellfun(@(x)x(keepIndices,:),data.subsegmentsList,'UniformOutput',false);
data.neighborIDs = data.neighborIDs(keepIndices,:);
interfaceLabels = interfaceLabels(keepIndices);
undecidedList = undecidedList(keepIndices);


end

