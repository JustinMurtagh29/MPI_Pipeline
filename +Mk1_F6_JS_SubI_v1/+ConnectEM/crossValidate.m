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
addpath('/gaba/u/sahilloo/repos/benedikt', '-end');
addpath(genpath('/gaba/u/sahilloo/repos/amotta/matlab/'))

%{
%% Configuration
config = 'ex144_08x2_mrNet';

nmlDir = Util.getGitReposOnPath();
nmlDir = fullfile( ...
    nmlDir{1}, 'data', 'tracings', ...
    'ex144-08x2-mrNet', 'connect-em');
%}

rootDir = '/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet/';
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

nmlDir = fullfile( param.saveFolder, 'tracings', 'connectEM','proofread');

folds = 5;
validSize = 200;
trainSizes = [200, 400, 600, 800, inf];

warning('Constructing FeatureMap de-novo');
fm = SynEM.getFeatureMap('paper');

info = Util.runInfo();
Util.showRunInfo(info);
%{
%% Loading data
config = loadConfig(config);
param = config.param;
%}
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

classifiers = cell(0);
results = cell(0, 2);
for curTrainSize = trainSizes
    for curFold = 1:folds
        rng(curFold);
        
        curRandIds = randperm(height(gt));
        curValIds = curRandIds(1:validSize);
        
        curTrainIds = min(validSize + [1, curTrainSize], height(gt));
        curTrainIds = curRandIds(curTrainIds(1):curTrainIds(2));
        
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
        curValLabels = repmat(gt.label(curValIds), 2, 1);
        
       [~, curScores] = predict(curClassifier, curValFeats);
        curScores = curScores(:, 1);
        
        classifiers{end + 1} = curClassifier; %#ok
        results(end + 1, :) = {curValLabels, curScores}; %#ok
    end
end

classifiers = reshape( ...
    classifiers, [folds, numel(trainSizes)]);
results = reshape( ...
    results, [folds, numel(trainSizes), 2]);

%% Plot aggregate results
clear cur*;

curFig = figure();
curFig.Color = 'white';

curAx = axes(curFig);
curAx.TickDir = 'out';
hold(curAx, 'on');

for curTrainSizeIdx = 1:numel(trainSizes)
    curTrainSize = trainSizes(curTrainSizeIdx);
    
    curLabels = results(:, curTrainSizeIdx, 1);
    curScores = results(:, curTrainSizeIdx, 2);
    
   [~, curIds] = cellfun( ...
        @sort, curScores, 'UniformOutput', false);
    curLabels = cellfun( ...
        @(l, s) l(s), curLabels, curIds, 'UniformOutput', false);
    curScores = cellfun( ...
        @(s) transpose(1:numel(s)), curScores, 'UniformOutput', false);
    
    curScores = cell2mat(curScores);
    curLabels = cell2mat(curLabels);
    
    [curPrec, curRec] = ...
        TypeEM.Classifier.buildPrecRec(curScores, curLabels);
    curLine = ...
        TypeEM.Classifier.plotPrecisionRecall(curAx, curPrec, curRec);
end

curLines = flip(curAx.Children);

xlim(curAx, [50, 100]);
ylim(curAx, [50, 100]);
axis(curAx, 'square');

curLegs = arrayfun( ...
    @(n) sprintf('%d training edges', n), ...
    trainSizes, 'UniformOutput', false);
curLegs = legend( ...
    curLines, curLegs, ...
    'Location', 'EastOutside');
curLegs.Box = 'off';

title(curAx, ...
    {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
savefig(fullfile(param.saveFolder,'connectEM','precrec_3boxes.fig'))
saveas(gcf,fullfile(param.saveFolder,'connectEM','precrec_3boxes.png'))
