load('E:\workspace\data\SynapseTrainingDataOverview.mat')
basefolder = 'E:\workspace\data\SynapseImageDataSet';
if ~exist(basefolder,'dir')
    mkdir(basefolder)
end
imSize = [61, 61];
% numClassInstances = 10;
augmentData = true;
imageFormat = 'jpg';
files = find([dataStruct.noVoxelLabel]);

%% make training set
trainingFolder = [basefolder filesep 'training'];
for s = files(1:end-1)
    fprintf('Processing file %s\n',dataStruct(s).file);
    m = matfile(dataStruct(s).file);
    data = m.data;
    raw = data.raw;
    tmp = m.voxelLabels;
    target = tmp(101:400,101:400,41:160);
    numClassInstances = floor(sum(target(:))*0.8);
    Codat.Caffe.makeDataset(raw, target, trainingFolder, imSize, numClassInstances, augmentData, imageFormat);
end

%% make test set
trainingFolder = [basefolder filesep 'test'];
for s = files(end)
    fprintf('Processing file %s\n',dataStruct(s).file);
    m = matfile(dataStruct(s).file);
    data = m.data;
    raw = data.raw;
    tmp = m.voxelLabels;
    target = tmp(101:400,101:400,41:160);
    numClassInstances = floor(sum(target(:))*0.8);
    Codat.Caffe.makeDataset(raw, target, trainingFolder, imSize, numClassInstances, augmentData, imageFormat);
end