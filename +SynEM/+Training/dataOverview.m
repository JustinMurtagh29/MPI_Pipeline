function info = dataOverview( dataFolder, areaThreshold )
%DATAOVERVIEW Generate data overview.
% INPUT dataFolder: Parent folder of synapse detection training data.
%       areaThresold: Optional parameter specifying an area threshold to
%       	use. All interfaces greater than the threshold are kept.
%       	(Default: 150)
% OUTPUT info: Table containing information about all training cubes in the
%           specified dataFolder.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

dataFolder = Util.addFilesep(dataFolder);
if ~exist('areaThreshold','var') || isempty(areaThreshold)
    areaThreshold = 150;
end

valCubes = { ...
    '2012-09-28_ex145_07x2_mag1_x1501-2000_y1501-2000_z1501-1700.mat', ...
    '2012-09-28_ex145_07x2_mag1_x1501-2000_y1501-2000_z2001-2200.mat', ...
    '2012-09-28_ex145_07x2_mag1_x2001-2500_y1001-1500_z1001-1200.mat', ...
    '2012-09-28_ex145_07x2_mag1_x2001-2500_y2001-2500_z1901-2100.mat', ...
    '2012-09-28_ex145_07x2_mag1_x2001-2500_y3001-3500_z801-1000.mat', ...
    '2012-09-28_ex145_07x2_mag1_x3001-3500_y3001-3500_z1401-1600.mat', ...
    '2012-09-28_ex145_07x2_mag1_x4001-4500_y2501-3000_z1001-1200.mat', ...
    '2012-09-28_ex145_07x2_mag1_x7001-7500_y5001-5500_z1501-1700.mat'};

unusedCubes = { ...
    '2012-09-28_ex145_07x2_mag1_x2001-2500_y2501-3000_z1001-1200.mat', ...
    '2012-09-28_ex145_07x2_mag1_x2001-2500_y2501-3000_z2001-2200.mat', ...
    '2012-09-28_ex145_07x2_mag1_x2001-2500_y4001-4500_z1001-1200.mat', ...
    '2012-09-28_ex145_07x2_mag1_x2001-2500_y5001-5500_z1001-1200.mat', ...
    '2012-09-28_ex145_07x2_mag1_x2001-2500_y501-1000_z1401-1600.mat', ...
    '2012-09-28_ex145_07x2_mag1_x4501-5000_y2501-3000_z401-600.mat', ...
    '2012-09-28_ex145_07x2_mag1_x4501-5000_y501-1000_z801-1000.mat'};

s = what(dataFolder);
Util.log('Inspecting %d files in %s.', length(s.mat), dataFolder);
info = struct;
for file = 1:length(s.mat)
    [~, metadata, interfaceLabels, ~, voxelLabels] = loadInterfaces( ...
        fullfile(dataFolder, s.mat{file}), areaThreshold);
    info(file).file = fullfile(dataFolder, s.mat{file});
    info(file).bboxBig = metadata.bboxBig;
    info(file).bboxSmall = metadata.bboxSmall;
    info(file).noInterfaces = length(interfaceLabels);
    info(file).noSynapses = sum(interfaceLabels < 3);
    info(file).noVoxelLabels = sum(voxelLabels(:));
    info(file).areaThreshold = areaThreshold;
    info(file).usage = 0;
    if any(strcmp(s.mat{file},valCubes))
        info(file).usage = 1;
    elseif any(strcmp(s.mat{file},unusedCubes))
        info(file).usage = 2;
    end
    fprintf('.'); %to see some progress
end
fprintf('\n');
Util.log('Finished training file inspection.');

info = struct2table(info);

%convert usage to categorical
info.usage = categorical(info.usage, 0:2, ...
    {'training','validation','unused'});

end

