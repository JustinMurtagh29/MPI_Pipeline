% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
% Modified by
%   Sahil Loomba <sahil.loomba@brain.mpg.de>
clear;
methodUsed = 'RUSBoost'; %'AdaBoostM1'; % 'LogitBoost';
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

% train with increasing training data sizes
trainFrac = [0.2, 0.4, 0.6, 0.8, 1];
trainSizes = floor(trainFrac.*height(gtTrain));
rng(1); % for reproducibility
curRandIds = randperm(height(gtTrain));

classifiers = cell(0);
results = cell(0, 2);

% compile test data
curTestFeats = [ ...
   [rawFeatsTest; ...
    fm.invertDirection(rawFeatsTest)], ...
   [classFeatsTest; ...
    fm.invertDirection(classFeatsTest)]];
curTestLabels = repmat(gtTest.label, 2, 1);
gtTestModified = repmat(gtTest,2,1); % because of inverted features

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

       [~, curScores] = predict(curClassifier, curTestFeats);
        curScores = curScores(:, 1);
        
        classifiers{end + 1} = curClassifier; %#ok
        results(end + 1, :) = {curTestLabels, curScores}; %#ok
end
%% Plot aggregate results
curFig = figure();
curFig.Color = 'white';

curAx = axes(curFig);
curAx.TickDir = 'out';
hold(curAx, 'on');

for curTrainSizeIdx = 1:numel(trainSizes)
    curTrainSize = trainSizes(curTrainSizeIdx);
    
    curLabels = results(curTrainSizeIdx, 1);
    curScores = results(curTrainSizeIdx, 2);
    
    curScores = cell2mat(curScores);
    curLabels = cell2mat(curLabels);
    
    [curPrec, curRec, threshVec] = ...
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

saveas(gcf,fullfile(param.saveFolder,'connectEM',['precrec_test_' methodUsed '.png']))

%% Building output
Util.log('Building output');
clear cur*;

classifier = struct;
classifier.gt = gt;
classifier.classifier = classifiers{1}; %choose trained on first set of edges
classifier.info = info;

Util.save(['/u/sahilloo/H2_3_v2_U1_SubI/connect-em/edgeClassifier/' datestr(clock,30) '.mat'],classifier);

%{
%% Inspect FPs
param.experimentName = 'H2_3_v2_U1_SubI_mr2e_wsmrnet';
gtTestModified = repmat(gtTest,2,1); % because of inverted features
outDir = fullfile(param.saveFolder,'connectEM');

scoreThr = -15.0287; % 75% recall, 93% precision
skel = H2_3_v2_SubI.ConnectEM.inspectFPs(param, gtTestModified, curScores, scoreThr);
skel.write(fullfile(outDir, sprintf('fps-%03f.nml',scoreThr)));
%}
