function dataStruct = generateDataOverview( dataFolder, areaThreshold )
%GENERATEDATAOVERVIEW Generate data overview.
% INPUT dataFolder: Parent folder of synapse detection training data.
%       areaThresold: Optional parameter specifying an area threshold to
%                     use. All interfaces greater than the threshold are
%                     kept. (Default: 150)
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('areaThreshold','var') || isempty(areaThreshold)
    areaThreshold = 150;
end

s = what(dataFolder);
dataStruct = struct;
for file = 1:length(s.mat)
    [~, metadata, interfaceLabels, ~, voxelLabels] = loadInterfaces([dataFolder, filesep, s.mat{file}],areaThreshold);
    dataStruct(file).file = [dataFolder, filesep, s.mat{file}];
    dataStruct(file).bboxBig = metadata.bboxBig;
    dataStruct(file).bboxSmall = metadata.bboxSmall;
    dataStruct(file).noInterfaces = length(interfaceLabels);
    dataStruct(file).noSynapses = sum(interfaceLabels < 3);
    dataStruct(file).noVoxelLabels = sum(voxelLabels(:));
    dataStruct(file).areaThreshold = areaThreshold;
end



end

