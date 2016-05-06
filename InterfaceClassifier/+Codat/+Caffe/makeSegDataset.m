%% load and set parameters
load('cortexTrainingDataParameter.mat')
% basefolder = 'E:\workspace\data\SegmentationDataSet';
basefolder = '/media/benedikt/DATA2/workspace/data/SegmentationDataSet';
if ~exist(basefolder,'dir')
    mkdir(basefolder)
end
imSize = [61, 61];
numClassInstances = 10000;
augmentData = true;
imageFormat = 'tif';
stackFiles = find(~exclude);

%% make training dataset
trainingFolder = [basefolder filesep 'training'];
for s = stackFiles(1:6)'
    fprintf('Processing stack file %d\n',s);
%     m = matfile(stacks(s).stackFile);
    m = matfile(strrep(stacks(s).stackFile,'E:\workspace\data\cortexTrainingData\targetKLEE\','/media/benedikt/DATA2/workspace/data/cortexTrainingData/targetKLEE/'));
    raw = m.raw;
    target = m.target;
    %set values in target
    targetTmp = zeros(size(target));
    targetTmp(target == 1) = 0; %intracellular
    targetTmp(target == -1) = 1; %extracellular
    targetTmp(target == 0) = -1; %discard
    Codat.Caffe.makeDataset(raw, targetTmp, trainingFolder, imSize, numClassInstances, augmentData, imageFormat);
end

%% make test set
validationFolder = [basefolder filesep 'test'];
for s = stackFiles(end-1:end)'
    fprintf('Processing stack file %d\n',s);
%     m = matfile(stacks(s).stackFile);
    m = matfile(strrep(stacks(s).stackFile,'E:\workspace\data\cortexTrainingData\targetKLEE\','/media/benedikt/DATA2/workspace/data/cortexTrainingData/targetKLEE/'));
    raw = m.raw;
    target = m.target;
    %set values in target
    targetTmp = zeros(size(target));
    targetTmp(target == 1) = 0; %intracellular
    targetTmp(target == -1) = 1; %extracellular
    targetTmp(target == 0) = -1; %discard
    Codat.Caffe.makeDataset(raw, targetTmp, validationFolder, imSize, numClassInstances, false, imageFormat);
end