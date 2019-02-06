% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% HACKHACKHACK
% NOTE(amotta): This is a huge mess. The training data is located in my
% repository, SynEM is from Benedikt's repository, and the SynEM classifier
% is loaded from the pipeline repo.
%   Let's get rid of the dependency on the pipeline repo by building the
% feature map de-novo. The `paper` version is identical to the one stored
% in the pipeline repository up to different `border` values in the
% `AverageFilter`. But this doesn't affect `fm.invertDirection`.
addpath('/gaba/u/amotta/code/benedikt', '-end');

%% Configuration
config = 'ex144_08x2_mrNet';

nmlDir = Util.getGitReposOnPath();
nmlDir = fullfile( ...
    nmlDir{1}, 'data', 'tracings', ...
    'ex144-08x2-mrNet', 'connect-em');

% Fraction of ground-truth to use for validation
validFrac = 0.1;

warning('Constructing FeatureMap de-novo');
fm = SynEM.getFeatureMap('paper');

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
config = loadConfig(config);
param = config.param;

%% Load ground truth
Util.log('Loading ground truth annotations');

curFiles = dir(fullfile(nmlDir, '*.nml'));
curFiles(cat(1, curFiles.isdir)) = [];
curFiles = fullfile(nmlDir, {curFiles.name});

gt = ConnectEM.loadGroundTruth(curFiles);

% Make sure that there are no multiplicities
assert(numel(gt.borderId) == numel(unique(gt.borderId)));

% Ignore borders without labels
gt(~gt.label, :) = [];

%% Load feature data for all borders
Util.log('Loading edge features');
clear cur*;

[rawFeats, rawBorderIds] = ...
    ConnectEM.loadFeaturesForBorderIds( ...
        param, 'InterfaceRawFeatures.mat', gt.borderId);
[classFeats, classBorderIds] = ...
    ConnectEM.loadFeaturesForBorderIds( ...
        param, 'InterfaceClassFeatures.mat', gt.borderId);
assert(isequal(rawBorderIds, classBorderIds));

[~, curRowIds] = ismember(rawBorderIds, gt.borderId);
gt = gt(curRowIds, :);

%% Separate into training and validation set
clear cur*;

rng(0);
curRandIds = randperm(height(gt));

curValidCount = round(validFrac * height(gt));
validIds = reshape(curRandIds(1:curValidCount), [], 1);
trainIds = reshape(curRandIds((1 + curValidCount):end), [], 1);

% Sanity check
assert(isempty(intersect(trainIds, validIds)));

%% Train classifier
Util.log('Training classifier');
clear cur*;

curFeats = [ ...
   [rawFeats(trainIds, :); ...
    fm.invertDirection(rawFeats(trainIds, :))], ...
   [classFeats(trainIds, :); ...
    fm.invertDirection(classFeats(trainIds, :))]];
curLabels = repmat(gt.label(trainIds), 2, 1);

classifier = fitensemble( ...
    curFeats, curLabels, ...
    'LogitBoost', 1500, 'tree', ...
    'LearnRate', 0.1, ...
    'NPrint', 100, ...
    'Prior', 'Empirical', ...
    'Type', 'Classification', ...
    'ClassNames', [+1, -1]);

%% Evaluation of classifier
Util.log('Evaluating classifier');
clear cur*;

curFeats = [ ...
   [rawFeats(validIds, :); ...
    fm.invertDirection(rawFeats(validIds, :))], ...
   [classFeats(validIds, :); ...
    fm.invertDirection(classFeats(validIds, :))]];

[~, validScores] = predict(classifier, curFeats);
validScores = reshape(validScores(:, 1), [], 2);

%% Building output
Util.log('Building output');
clear cur*;

out = struct;
out.gt = gt;
out.trainIds = trainIds;
out.validIds = validIds;
out.validScores = validScores;
out.classifier = classifier;
out.info = info;
