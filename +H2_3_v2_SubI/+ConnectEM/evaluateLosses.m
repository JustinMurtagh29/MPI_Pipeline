% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
% Modified by
%   Sahil Loomba <sahil.loomba@brain.mpg.de>
methodUsed = 'LogitBoost';
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

rootDir = '/tmpscratch/sahilloo/data/H2_3_v2_U1_SubI/pipelineRun_mr2e_wsmrnet/';
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

nmlDir = fullfile( param.saveFolder, 'tracings', 'connectEM','proofread');

warning('Constructing FeatureMap de-novo');
fm = SynEM.getFeatureMap('paper');

info = Util.runInfo();
Util.showRunInfo(info);

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
Util.log('Done loading features and labels...')

%% divide into training and test set
rng(0);
testFrac = 0.2;
testSize = ceil(testFrac*height(gt));
curRandIds = randperm(height(gt));
curTestIds = curRandIds(1:testSize);
curTrainIds = curRandIds(testSize+1:end);
% test features
gtTest = gt(curTestIds,:);
rawFeatsTest = rawFeats(curTestIds,:);
classFeatsTest = classFeats(curTestIds,:);
% train features
gtTrain = gt(curTrainIds,:);
rawFeatsTrain = rawFeats(curTrainIds,:);
classFeatsTrain = classFeats(curTrainIds,:);

clear cur*

% compile test data
curTestFeats = [ ...
   [rawFeatsTest; ...
    fm.invertDirection(rawFeatsTest)], ...
   [classFeatsTest; ...
    fm.invertDirection(classFeatsTest)]];
curTestLabels = repmat(gtTest.label, 2, 1);

% train with increasing training data sizes
trainFrac = [0.2, 0.4, 0.6, 0.8, 1];
trainSizes = floor(trainFrac.*height(gtTrain));
rng(1); % for reproducibility
curRandIds = randperm(height(gtTrain));

classifiers = cell(0);

curFig = figure();
curFig.Color = 'white';
curAx = axes(curFig);
hold(curAx,'on');
for curTrainSize = trainSizes
        curTrainIds = curRandIds(1:curTrainSize);
        
        curTrainFeats = [ ...
           [rawFeatsTrain(curTrainIds, :); ...
            fm.invertDirection(rawFeatsTrain(curTrainIds, :))], ...
           [classFeatsTrain(curTrainIds, :); ...
            fm.invertDirection(classFeatsTrain(curTrainIds, :))]];
        curTrainLabels = repmat(gtTrain.label(curTrainIds), 2, 1);
        
        curClassifier = fitensemble( ...
            curTrainFeats, curTrainLabels, ...
            methodUsed, 1500, 'tree', ...
            'LearnRate', 0.1, ...
            'NPrint', 100, ...
            'Prior', 'Empirical', ...
            'Type', 'Classification', ...
            'ClassNames', [+1, -1]);
        
        classifiers{end + 1} = curClassifier; %#ok
        rsLoss = resubLoss(curClassifier,'Mode','Cumulative');
        plot(rsLoss);
        ylabel('Resubstitution Loss');
        xlabel('Number of Learning Cycles');
end
curLines = flip(curAx.Children);
curLegs = arrayfun( ...
    @(n) sprintf('%d training edges', n), ...
    trainSizes, 'UniformOutput', false);

curLegs = legend( ...
    curLines, curLegs, ...
    'Location', 'EastOutside');
curLegs.Box = 'off';
saveas(gcf,fullfile(param.saveFolder,'connectEM','validationResubLoss.png'))
close all

%% with holdout
classifiers = cell(0);

curFig = figure();
curFig.Color = 'white';
curAx = axes(curFig);
hold(curAx,'on');
for curTrainSize = trainSizes
        curTrainIds = curRandIds(1:curTrainSize);

        curTrainFeats = [ ...
           [rawFeatsTrain(curTrainIds, :); ...
            fm.invertDirection(rawFeatsTrain(curTrainIds, :))], ...
           [classFeatsTrain(curTrainIds, :); ...
            fm.invertDirection(classFeatsTrain(curTrainIds, :))]];
        curTrainLabels = repmat(gtTrain.label(curTrainIds), 2, 1);

        curClassifier = fitensemble( ...
            curTrainFeats, curTrainLabels, ...
            methodUsed, 1500, 'tree', ...
            'LearnRate', 0.1, ...
            'NPrint', 100, ...
            'Prior', 'Empirical', ...
            'Type', 'Classification', ...
            'ClassNames', [+1, -1],...
            'Holdout',0.2);

        classifiers{end + 1} = curClassifier; %#ok
        genError = kfoldLoss(curClassifier,'Mode','Cumulative');
        plot(genError);
        ylabel('Generalization Error');'RUSBoost'
        xlabel('Number of Learning Cycles');
end
curLines = flip(curAx.Children);
curLegs = arrayfun( ...
    @(n) sprintf('%d training edges', n), ...
    trainSizes, 'UniformOutput', false);

curLegs = legend( ...
    curLines, curLegs, ...
    'Location', 'EastOutside');
curLegs.Box = 'off';
saveas(gcf,fullfile(param.saveFolder,'connectEM',['validationGenError' methodUsed '.png'))
