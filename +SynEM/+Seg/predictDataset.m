function [p, job] = predictDataset( p, fm, classifier, outputFile, ...
    cluster, cubeIdx, featureFile, interfaceFile, normalizeRaw, ...
    saveClassAndFM )
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
%       featureFile: (Optional) logical or string
%           Flag indicating whether to save features for each local
%           segmentation cube. If a string is specified this string will be
%           used as the name of the feature file in the local cubes. The
%           default name is 'InterfaceFeatures.mat'.
%           (Default: false)
%       interfaceFile: (Optional) logical or string
%           Flag indicating whether to save interfaces for each local
%           segmentation cube. If a string is specified this string will be
%           used as the name of the interface file in the local cubes. The
%           default name is 'Interfaces.mat'.
%           (Default: false)
%       normalizeRaw: (Optional) logical
%           Flag to indicate that the raw data is normalized to 122 mean
%           and 22 standard deviation by using p.norm.func to normalize to
%           mean 0 and 1 std first.
%           (Default: no raw data normalization).
%       saveClassAndFM: (Optional) string
%           Name of the file with classifier and feature map in the
%           segmentation main folder.
%           (default: 'SynapseClassifier.mat')
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

info = Util.runInfo(true);

if ~exist('outputFile','var') || isempty(outputFile)
    outputFile = 'synapseScores';
end
[~,  outputFile] = fileparts(outputFile);
outputFile = [outputFile, '.mat'];

% parse cube idx input
if ~exist('cubeIdx','var') || isempty(cubeIdx)
    cubeIdx = 1:numel(p.local);
elseif islogical(cubeIdx)
    cubeIdx = find(cubeIdx);
end
Util.log('Running synapse detection for %d local cubes.', length(cubeIdx));

% parse input for feature file
if ~exist('featureFile','var') || isempty(featureFile)
    featureFile = [];
elseif ischar(featureFile)
    [~, name] = fileparts(featureFile);
    featureFile = [name '.mat'];
elseif featureFile
    featureFile = 'InterfaceFeatures.mat';
else
    featureFile = [];
end
if ischar(featureFile)
    Util.log('Storing features at %s.', featureFile);
else
    Util.log('Features are not being stored.');
end

% parse input for interface file
if ~exist('interfaceFile','var') || isempty(interfaceFile)
    interfaceFile = [];
elseif ischar(interfaceFile)
    [~, name] = fileparts(interfaceFile);
    interfaceFile = [name '.mat'];
elseif interfaceFile
    interfaceFile = 'Interfaces.mat';
else
    interfaceFile = [];
end
if ischar(featureFile)
    Util.log('Storing interfaces at %s.', featureFile);
else
    Util.log('Interfaces are not being stored.');
end

if ~exist('normalizeRaw','var') || isempty(normalizeRaw)
    normalizeRaw = false;
end
Util.log('Raw data normalization is set to %d.', normalizeRaw);

cubeIdx = cubeIdx(:)';

% get cluster object
if ~exist('cluster','var') || isempty(cluster)
    try
        cluster = getCluster('cpu');
    catch
        cluster = parcluster();
    end
end
try
    info.param.cluster = cluster.IndependentSubmitFcn(2:end);
catch
    info.param.cluster = [];
end

% save classification data
if ~exist('saveClassAndFM', 'var') || isempty(saveClassAndFM)
    saveClassAndFM = 'SynapseClassifier.mat';
else
    [~, name] = fileparts(saveClassAndFM);
    saveClassAndFM = [name '.mat'];
end
p.synEM = [p.saveFolder saveClassAndFM];
Util.log('Saving classifier and feature map to %s.', p.synEM);
Util.save([p.synEM], classifier, fm, info);
    
% submit to cluster
Util.log('Submitting jobs to cluster.');
for i = cubeIdx
    p.local(i).synapseFile = [p.local(i).saveFolder, outputFile];
    if ~isempty(featureFile)
        p.local(i).interfaceFeatureFile = ...
            [p.local(i).saveFolder, featureFile];
    end
    if ~isempty(interfaceFile)
        p.local(i).interfaceFile = ...
            [p.local(i).saveFolder, interfaceFile];
    end
end

inputCell = arrayfun(@(i){[p.saveFolder 'allParameter.mat'], i, ...
    saveClassAndFM, outputFile, featureFile, interfaceFile, ...
    normalizeRaw}, cubeIdx, ...
    'UniformOutput', false);

job = SynEM.Util.startJob(cluster, @jobWrapper, inputCell, 0, ...
    'SynapseDetection_Job%d');
Util.log('Interface classification submitted as job %d.', job.Id);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function jobWrapper(allParamFile, i, classFile, outputFile, ...
    featureFile, interfaceFile, normalizeRaw)

tic
%load parameter file
m = load(allParamFile, 'p');
p = m.p;

%check for output
outputFile = [p.local(i).saveFolder outputFile];
if exist(outputFile, 'file')
    % error('Output file %s already exists.', outputFile )
end

%load fm and classifier
m = load([p.saveFolder classFile], 'fm', 'classifier');
fm = m.fm;
classifier = m.classifier;

%actual prediction
[scores, X, interfaces] = SynEM.Seg.predictCube(p, i, fm, classifier, ...
    normalizeRaw);
scores = scores(:,1); %for default matlab classifiers
if strcmp(fm.mode, 'direction')
    %both direction of one interface in a row
    scores = reshape(scores, [], 2);
end

%save results
Util.log('Saving scores to %s.', outputFile);
runTime = toc;
Util.save(outputFile, scores, runTime);

if ~isempty(featureFile)
    if strcmp(fm.mode, 'direction')
        X = X(1:end/2,:); %only save first direction
    end
    outputFile = [p.local(i).saveFolder featureFile];
    Util.log('Saving features to %s.', outputFile);
    Util.save(outputFile, X);
end
if ~isempty(interfaceFile)
    outputFile = [p.local(i).saveFolder interfaceFile];
    Util.log('Saving interfaces to %s.', outputFile);
    Util.save(outputFile, interfaces);
end

end
