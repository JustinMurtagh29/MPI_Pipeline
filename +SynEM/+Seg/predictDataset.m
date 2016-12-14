function [p, job] = predictDataset( p, fm, classifier, outputFile, ...
    cluster, cubeIdx )
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
%           (Default: getCluster('cpu');
%       cubeIdx: (optional) [1xN] int
%           Linear or logical indices of the cubes in p.local for which the
%           prediction is done.
%           (Default: 1:numel(p.local))
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
if iscolumn(cubeIdx); cubeIdx = cubeIdx'; end

%save classification data
p.synEM = [p.saveFolder 'synapseClassifier.mat'];
Util.save(p.synEM, classifier, fm);
[~,  outputFile] = fileparts(outputFile);
outputFile = [outputFile, '.mat'];
for i = cubeIdx
    p.local(i).synapseFile = [p.local(i).saveFolder, outputFile];
end

inputCell = arrayfun(@(i){p, i, fm, classifier, outputFile}, ...
    cubeIdx, 'UniformOutput', false);

if ~exist('cluster','var') || isempty(cluster)
    cluster = getCluster('cpu');
end
job = SynEM.Util.startJob(cluster, @jobWrapper, inputCell, 0, ...
    'SynapseDetection');
end

function jobWrapper(p, i, fm, classifier, outputFile)

scores = SynEM.Seg.predictCube(p, i, fm, classifier);
scores = scores(:,1); %for default matlab classifiers
if strcmp(fm.mode,'direction')
    %both direction of one interface in a row
    scores = reshape(scores,[],2); 
end
outputFile = [p.local(i).saveFolder outputFile];
fprintf(['[%s] SynEM.Seg.predictCube - Saving output to ' ...
    'cube %s.\n'], datestr(now), outputFile);
Util.save(outputFile,scores);

end
