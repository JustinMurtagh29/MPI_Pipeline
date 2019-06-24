% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
% Modified by
%   Sahil Loomba <sahil.loomba@brain.mpg.de>
clear;
methodUsed = 'LogitBoost'; %'AdaBoostM1'; % 'LogitBoost';
addpath(genpath('/gaba/u/sahilloo/repos/amotta/matlab/'))

rootDir = '/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet/';
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;
param.experimentName = 'Mk1_F6_JS_SubI_v1_mrnet_wsmrnet';
import Mk1_F6_JS_SubI_v1.TypeEM.*

info = Util.runInfo();

% load training data
featureSetName = 'segmentAgglomerate';
nmlDir = fullfile(param.saveFolder, ...
     'tracings', 'typeEM');
nmlFiles = fullfile(nmlDir, 'proofread', ...
     {'box-1.nml','box-2.nml', 'box-3.nml', ...
     'box-4.nml','box-5.nml','box-6.nml',...
     'box-7.nml','box-8.nml','box-9.nml',...
     'box-10.nml','box-11.nml','box-12.nml',...
     'box-13.nml'});

rng(0);
idxTrain = [1,2,3,5,6,7,8,9,10,11,12,13];
idxTest = 4;
% load features
gt = TypeEM.GroundTruth.loadSet( ...
        param, featureSetName, nmlFiles(idxTrain));
gt = TypeEM.GroundTruth.loadFeatures(param, featureSetName, gt);
% load test set
gtTest = TypeEM.GroundTruth.loadSet( ...
        param, featureSetName, nmlFiles(idxTest));

curRandIds = randperm(size(gt.label,1));
trainIds = curRandIds;

% train with increasing training data sizes
trainFrac = [0.2, 0.4, 0.6, 0.8, 1];
trainSizes = floor(trainFrac.*length(trainIds));
curRandIds = randperm(length(trainIds));

classifiers = cell(0);
results = cell(0, 2);

for curTrainSize = trainSizes
    curTrainIds = curRandIds(1:curTrainSize);

    gtTrain = gt;
    gtTrain.segId = gt.segId(curTrainIds,:);
    gtTrain.label = gt.label(curTrainIds,:);
    gtTrain.featMat = gt.featMat(curTrainIds,:);

    curClassifier.classes = gt.class;
    curClassifier.classifiers = arrayfun( ...
        @(c) buildForClass(gtTrain, c, methodUsed), ...
        gtTrain.class, 'UniformOutput', false);
    curClassifier.featureSetName = featureSetName;

    % apply classifier to test data
    [precRec, fig, curGtTest] = TypeEM.Classifier.evaluate(param, curClassifier, gtTest);
    title([methodUsed ' with trainSize:' num2str(curTrainSize)])
    saveas(gcf,fullfile(param.saveFolder,'typeEM',['precrec_box_' methodUsed '_tsize_' num2str(curTrainSize) '.png']))
    close all

    % build platt parameters
    trainPlattFunc = @(curIdx) trainPlattForClass(curGtTest, curIdx);
    [aVec, bVec] = arrayfun(trainPlattFunc, 1:numel(curClassifier.classes));

    % build output
    platt = struct( ...
        'a', num2cell(aVec), ...
        'b', num2cell(bVec));
    % convert to probs
    curClassifier.plattParams = platt;
    probs = TypeEM.Classifier.applyPlatt(curClassifier, curGtTest.scores);
    curGtTest.probs = probs;

    classifiers{end + 1} = curClassifier;
    results(end + 1, :) = {precRec, probs};
end 
%{
% Look at false positives
curCount = 40;
className = 'glia';
skels = Debug.inspectFPs(param, curCount, className, curGtTest);
skels.write(fullfile(param.saveFolder,'typeEM', sprintf('fps-%s.nml',className)));
className = 'axon';
skels = Debug.inspectFPs(param, curCount, className, curGtTest);
skels.write(fullfile(param.saveFolder,'typeEM', sprintf('fps-%s.nml',className)));
className = 'dendrite';
skels = Debug.inspectFPs(param, curCount, className, curGtTest);
skels.write(fullfile(param.saveFolder,'typeEM', sprintf('fps-%s.nml',className)));

% Look at true positives
curCount = 100;
className = 'glia';
skels = Debug.inspectTPs(param, curCount, className, curGtTest);
skels.write(fullfile(param.saveFolder,'typeEM', sprintf('tps-%s.nml',className)));
className = 'axon';
skels = Debug.inspectTPs(param, curCount, className, curGtTest);
skels.write(fullfile(param.saveFolder,'typeEM', sprintf('tps-%s.nml',className)));
className = 'dendrite';
skels = Debug.inspectTPs(param, curCount, className, curGtTest);
skels.write(fullfile(param.saveFolder,'typeEM', sprintf('tps-%s.nml',className)));

% label statistics
Debug.labelDistributions

Util.save(['/u/sahilloo/Mk1_F6_JS_SubI_v1/type-em/typeClassifier/' datestr(clock,30) '.mat'],classifiers,gt,info);
%}
%{
%% Building output
Util.log('Building output');
classifier = classifiers{end}; %choose trained on all data
Util.save(['/u/sahilloo/Mk1_F6_JS_SubI_v1/type-em/typeClassifier/' datestr(clock,30) '.mat'],classifier,gt,info);
%}

function classifier = buildForClass(gt, class, methodUsed)
    % classifier = buildForClass(data, className)
    %   Builds a one-versus-all (OVA) classifier for a class

    classIdx = find(gt.class == class);
    mask = gt.label(:, classIdx) ~= 0;

    % extract training data
    trainFeat = gt.featMat(mask, :);
    trainClass = gt.label(mask, classIdx) > 0;

    % train classifier
    classifier = fitensemble( ...
        trainFeat, trainClass, ...
        methodUsed, 1500, 'tree', ...
        'LearnRate', 0.1, ...
        'NPrint', 100, ...
        'Prior', 'Empirical', ...
        'Type', 'Classification', ...
        'ClassNames', [true, false]);
end

function [a, b] = trainPlattForClass(gt, classIdx)
    mask = gt.label(:, classIdx) ~= 0;
    scores = gt.scores(mask, classIdx);
    labels = double(gt.label(mask, classIdx) > 0);

    % optimize Platt parameter
    [a, b] = TypeEM.Classifier.learnPlattParams(scores, labels);
end

