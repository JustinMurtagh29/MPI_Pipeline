function makeTrainingData( dataOverview, targetSize, border, targetDir )
%MAKETRAININGDATA Create training data and output stack files.
% INPUT dataOverview: Synapse training data parameter struct.
%       targetSize: Size of target cubes.
%       border: Size of border around target cubes in raw.
%       targetDir: Directory for output files without filesep at the end.
%                  (Default: Current directory.)
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('targetSize','var') || isempty(targetSize)
    targetSize = [100, 100, 60];
end

if ~exist('border','var') || isempty(border)
    border = [100, 100, 40];
end

if ~exist('targetDir','var') || isempty(targetDir)
    targetDir = pwd;
end

stacks = struct;
count = 1;
for cube = 1:length(dataOverview)
    if dataOverview(cube).noVoxelLabels > 0
        m = matfile(dataOverview(cube).file);
        data = m.data;
        target = m.voxelLabels;
        target = target(101:400,101:400,41:160);
        [x_train,y_train] = Codat.CNN.Misc.tileTrainingCubes(data.raw(1:500,1:500,1:200), target, targetSize, border);
        for i = 1:length(y_train)
            raw = x_train{i};
            target = y_train{i};
            if sum(target(:) == 1) > 100
                save([targetDir filesep num2str(count), '.mat'],'raw','target');
                stacks(count).stackFile = [targetDir filesep num2str(count), '.mat'];
                stacks(count).synapticVoxels = sum(target(:) == 1);
                count = count + 1;
            end
        end
    end
end
save([targetDir filesep 'synapseVoxelTrainingDataParameter.mat'],'stacks');


end

