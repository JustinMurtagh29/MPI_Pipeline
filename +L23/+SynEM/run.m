% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/tmpscratch/amotta/l23/2018-10-09-mrnet-pipeline-run';

% NOTE(amotta): This is a version of SynEM that only uses hand-designed
% feature and that was trained on the extended training set. According to
% Benedikt's PhD thesis, this is the configuration that performs best on
% the L2/3 dataset.
classifierFile = 'SynapseClassifier_synem_tr2.mat';
scoreFile = strrep(classifierFile, 'Classifier', 'Scores');
paramFile = strrep(classifierFile, 'Classifier', 'Parameters');

info = Util.runInfo();

%% Load data
param = fullfile(rootDir, 'allParameter.mat');
param = load(param);
param = param.p;

classifierFile = fullfile(rootDir, classifierFile);
paramFile = fullfile(rootDir, paramFile);

curData = load(classifierFile, 'classifier', 'fm');
classifier = curData.classifier;
fm = curData.fm;

%% Run synapse detection
clear cur*;

cluster = Cluster.config( ...
    'memory', 24, 'time', '100:00:00');
[synParam, job] = SynEM.Seg.predictDataset( ...
    param, fm, classifier, scoreFile, cluster, ...
    [], true, true, true, [], false);
wait(job);

%% Save synapse parameters
curOut = struct;
curOut.info = info;
curOut.p = synParam;

Util.saveStruct(paramFile, curOut);
Util.protect(paramFile);
