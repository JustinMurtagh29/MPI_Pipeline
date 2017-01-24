function [p, job] = predictDataset( p, fm, classifier, outputFile, ...
    cluster, cubeIdx, saveFeatures, saveInterfaces )
%PREDICTDATASET SynEM prediction for a dataset.
% INPUT p: struct
%           SegEM segmentation parameter struct.
%       fm: SynEM.FeatureMap
%           The feature map for prediction.
%       classifier: object or string
%           Classifier object (e.g. SynEM.Classifier) or path to classifier
%           object mat-file containing a 'classifier' variable.
%       outputFile: (optional) string
%           Name of the output file that is stored in each segmentation
%           cube. Synapse scores are saved in the output file as a [Nx1] or
%           [Nx2] float array, where N equals the length of the border list
%           of the corresponding cube with borders of
%           size > fm.areaThreshold. If fm.mode is 'direction' and the rows
%           in scores correspond to the two interface directions (with the
%           first direction being equal to the direction of the current
%           entry in edges).
%           (Default: 'synapseScores.mat')
%       cluster: (optional) parallel.cluster object
%           Cluster object to submit jobs.
%           (Default: parcluster())
%       cubeIdx: (optional) [1xN] int
%           Linear or logical indices of the cubes in p.local for which the
%           prediction is done.
%           (Default: 1:numel(p.local))
%       saveFeatures: (Optional) logical
%           Flag indicating whether to save features for each local
%           segmentation cube.
%           (Default: false)
%       saveInterfaces: (Optional) logical
%           Flag indicating whether to save interfaces for each local
%           segmentation cube.
%           (Default: false)
% OUTPUT p: struct
%           Modified segmentaton parameter struct. Classifier and feature
%           maps are stored at 'p.synEM' and each local cube now contains a
%           'synapseFile' path with the save location of the synapse
%           scores.
%        job: parallel.job object.
%           Job array of the prediction jobs.
%
% NOTE The classifier and feature maps are saved to
%      [p.saveFolder 'synapseClassifier.mat'].
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('outputFile','var') || isempty(outputFile)
    outputFile = 'synapseScores';
end
if ~exist('cubeIdx','var') || isempty(cubeIdx)
    cubeIdx = 1:numel(p.local);
elseif islogical(cubeIdx)
    cubeIdx = find(cubeIdx);
end
if ~exist('saveFeatures','var') || isempty(saveFeatures)
    saveFeatures = false;
end
if ~exist('saveInterfaces','var') || isempty(saveInterfaces)
    saveInterfaces = false;
end
if iscolumn(cubeIdx); cubeIdx = cubeIdx'; end

%save classification data
p.synEM = [p.saveFolder 'SynapseClassifier.mat'];
Util.save([p.synEM], classifier, fm);
[~,  outputFile] = fileparts(outputFile);
outputFile = [outputFile, '.mat'];
for i = cubeIdx
    p.local(i).synapseFile = [p.local(i).saveFolder, outputFile];
    if saveFeatures
        p.local(i).interfaceFeatureFile = ...
            [p.local(i).saveFolder, 'InterfaceFeatures.mat'];
    end
    if saveInterfaces
        p.local(i).interfaceFile = ...
            [p.local(i).saveFolder, 'Interfaces.mat'];
    end
end

inputCell = arrayfun(@(i){[p.saveFolder 'allParameter.mat'], i, ...
    outputFile, saveFeatures, saveInterfaces}, cubeIdx, ...
    'UniformOutput', false);

if ~exist('cluster','var') || isempty(cluster)
    try
        cluster = getCluster('cpu');
    catch
        cluster = parcluster();
    end
end
job = SynEM.Util.startJob(cluster, @jobWrapper, inputCell, 0, ...
    'SynapseDetection');
end

function jobWrapper(allParamFile, i, outputFile, ...
    saveFeatures, saveInterfaces)

%load parameter file
m = load(allParamFile, 'p');
p = m.p;

%load fm and classifier
m = load([p.saveFolder 'SynapseClassifier.mat'], 'fm', 'classifier');
fm = m.fm;
classifier = m.classifier;

%actual prediction
[scores, X, interfaces] = SynEM.Seg.predictCube(p, i, fm, classifier);
scores = scores(:,1); %for default matlab classifiers
if strcmp(fm.mode, 'direction')
    %both direction of one interface in a row
    scores = reshape(scores, [], 2);
end

%save results
outputFile = [p.local(i).saveFolder outputFile];
fprintf('[%s] SynEM.Seg.predictCube - Saving output to %s.\n', ...
    datestr(now), outputFile);
Util.save(outputFile, scores);

if saveFeatures
    if strcmp(fm.mode, 'direction')
        X = X(1:end/2,:); %only save first direction
    end
    outputFile = [p.local(i).saveFolder 'InterfaceFeatures.mat'];
    fprintf('[%s] SynEM.Seg.predictCube - Saving features to %s.\n', ...
        datestr(now), outputFile);
    Util.save(outputFile, X);
end
if saveInterfaces
    outputFile = [p.local(i).saveFolder 'Interfaces.mat'];
    fprintf('[%s] SynEM.Seg.predictCube - Saving interfaces to %s.\n', ...
        datestr(now), outputFile);
    Util.save(outputFile, interfaces);
end

end
