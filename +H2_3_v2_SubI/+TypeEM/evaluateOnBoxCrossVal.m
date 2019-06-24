% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
% Modified by
%   Sahil Loomba <sahil.loomba@brain.mpg.de>
clear;
methodUsed = 'LogitBoost'; %'AdaBoostM1'; % 'LogitBoost';
addpath(genpath('/gaba/u/sahilloo/repos/amotta/matlab/'))

rootDir = '/tmpscratch/sahilloo/data/H2_3_v2_U1_SubI/pipelineRun_mr2e_wsmrnet/';
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;
param.experimentName = 'H2_3_v2_U1_SubI_mr2e_wsmrnet';
import H2_3_v2_SubI.TypeEM.*

info = Util.runInfo();

% load training data
featureSetName = 'segmentAgglomerate';
nmlDir = fullfile(param.saveFolder, ...
     'tracings', 'typeEM');
nmlFiles = fullfile(nmlDir, 'proofread', ...
     {'box-1.nml','box-2.nml', 'box-3.nml', ...
     'box-4.nml','box-5.nml','box-6.nml',...
     'box-7.nml','box-8.nml','box-9.nml'});

% do cross validation
classifiers = cell(0);
results = cell(0);

rng(0);
ids = 1:numel(nmlFiles);

for curBoxId=1:numel(nmlFiles)
idxTest = curBoxId;
idxTrain = setdiff(ids,curBoxId);

% load features
gtTrain = TypeEM.GroundTruth.loadSet( ...
        param, featureSetName, nmlFiles(idxTrain));
gtTrain = TypeEM.GroundTruth.loadFeatures(param, featureSetName, gtTrain);
% load test set
gtTest = TypeEM.GroundTruth.loadSet( ...
        param, featureSetName, nmlFiles(idxTest));

curClassifier.classes = gtTrain.class;
curClassifier.classifiers = arrayfun( ...
    @(c) buildForClass(gtTrain, c, methodUsed), ...
    gtTrain.class, 'UniformOutput', false);
curClassifier.featureSetName = featureSetName;

% apply classifier to test data
[precRec, fig, curGtTest] = TypeEM.Classifier.evaluate(param, curClassifier, gtTest);

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

saveas(gcf,fullfile(param.saveFolder,'typeEM',['precrec_box_crossVal' methodUsed '_' num2str(curBoxId) '.png']))
close all
end

% get train box statistics
numGlia = [];
numGliaN = [];
numAxon = [];
numAxonN = []
numDend = [];
numDendN = []
for i=1:numel(classifiers)
    thisG = classifiers{i}.classifiers{1};
    thisA = classifiers{i}.classifiers{2};
    thisD = classifiers{i}.classifiers{3};

    numGlia = vertcat(numGlia, sum(thisG.Y));
    numGliaN = vertcat(numGliaN, sum(~thisG.Y));
    numAxon = vertcat(numAxon, sum(thisA.Y));
    numAxonN = vertcat(numAxonN, sum(~thisA.Y));
    numDend = vertcat(numDend, sum(thisD.Y));
    numDendN = vertcat(numDendN, sum(~thisD.Y));
end


Util.save(fullfile(param.saveFolder,'typeEM','precrec_box_crossVal.mat'),numGlia, numGliaN, numAxon, numAxonN, numDend, numDendN)
figure
hold on
bar(cat(2,numGlia, numAxon, numDend))
bar(cat(2,numGliaN, numAxonN, numDendN),'FaceColor','none','EdgeColor','g')
legend({'+ glia','+ axon','+ dendrite','- glia','- axon','- dendrite' })
xlabel('Test box Id')
ylabel('Numder of labels')
saveas(gcf,fullfile(param.saveFolder,'typeEM','crossVal_trainBox_statistics.png'))
close all



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

