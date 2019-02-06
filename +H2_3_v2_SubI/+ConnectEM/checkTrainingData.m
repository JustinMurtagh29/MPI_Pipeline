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

folds = 3;

warning('Constructing FeatureMap de-novo');
fm = SynEM.getFeatureMap('paper');

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
config = loadConfig(config);
param = config.param;

%% Load ground truth
nmlFiles = dir(fullfile(nmlDir, '*.nml'));
nmlFiles(cat(1, nmlFiles.isdir)) = [];
nmlFiles = fullfile(nmlDir, {nmlFiles.name});

gt = ConnectEM.loadGroundTruth(nmlFiles);

% Make sure that there are no multiplicities
assert(numel(gt.borderId) == numel(unique(gt.borderId)));

% Ignore borders without labels
gt(~gt.label, :) = [];

%% Load feature data for all borders
[rawFeats, rawBorderIds] = ...
    ConnectEM.loadFeaturesForBorderIds( ...
        param, 'InterfaceRawFeatures.mat', gt.borderId);
[classFeats, classBorderIds] = ...
    ConnectEM.loadFeaturesForBorderIds( ...
        param, 'InterfaceClassFeatures.mat', gt.borderId);
assert(isequal(rawBorderIds, classBorderIds));

[~, curRowIds] = ismember(rawBorderIds, gt.borderId);
gt = gt(curRowIds, :);

%% Run cross-validation
clear cur*;

out = struct;
out.classifiers = cell(folds, 1);
out.validResults = cell(folds, 2);

rng(0);
curRandIds = randperm(height(gt));
curValidSize = ceil(height(gt) / folds);

for curFold = 1:folds
        curValIds = (curFold - [1, 0]) * curValidSize + [1, 0];
        curValIds = curRandIds(curValIds(1):min(curValIds(2), height(gt)));
        curTrainIds = setdiff(curRandIds, curValIds, 'stable');
        
        % Sanity check
        assert(isempty(intersect(curTrainIds, curValIds)));
        
        curTrainFeats = [ ...
           [rawFeats(curTrainIds, :); ...
            fm.invertDirection(rawFeats(curTrainIds, :))], ...
           [classFeats(curTrainIds, :); ...
            fm.invertDirection(classFeats(curTrainIds, :))]];
        curTrainLabels = repmat(gt.label(curTrainIds), 2, 1);
        
        curClassifier = fitensemble( ...
            curTrainFeats, curTrainLabels, ...
            'LogitBoost', 1500, 'tree', ...
            'LearnRate', 0.1, ...
            'NPrint', 100, ...
            'Prior', 'Empirical', ...
            'Type', 'Classification', ...
            'ClassNames', [+1, -1]);
        
        curValFeats = [ ...
           [rawFeats(curValIds, :); ...
            fm.invertDirection(rawFeats(curValIds, :))], ...
           [classFeats(curValIds, :); ...
            fm.invertDirection(classFeats(curValIds, :))]];
        
       [~, curScores] = predict(curClassifier, curValFeats);
        curScores = reshape(curScores(:, 1), [], 2);
        
        out.classifiers{curFold} = curClassifier;
        out.validResults(curFold, :) = {curValIds(:), curScores};
end

out.gt = gt;
out.info = info;

%% Export top 50 FP-edges per fold
clear cur*;

for curFold = 1:folds
    curGt = out.gt(out.validResults{curFold, 1}, :);
    curGt.score = out.validResults{curFold, 2};
    
    % Look at top false-positive candidates
    curGt = curGt(curGt.label < 0, :);
    curGt.score = max(curGt.score, [], 2);
    curGt = sortrows(curGt, 'score', 'descend');
end
