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

% load training and test data
featureSetName = 'segmentAgglomerate';
nmlDir = fullfile(param.saveFolder, ...
     'tracings', 'typeEM');
nmlFiles = fullfile(nmlDir, 'proofread', ...
     {'box-1.nml', 'box-2.nml', 'box-3.nml', ...
     'box-4.nml','box-5.nml','box-6.nml',...
     'box-7.nml','box-8.nml', 'box-9.nml'});
% load features
gt = TypeEM.GroundTruth.loadSet( ...
        param, featureSetName, nmlFiles);
gt = TypeEM.GroundTruth.loadFeatures(param, featureSetName, gt);

%% divide into training and test set
rng(0);
testFrac = 0.2;
testSize = ceil(testFrac*size(gt.label,1));
curRandIds = randperm(size(gt.label,1));
testIds = curRandIds(1:testSize);
trainIds = curRandIds(testSize+1:end);
gtTest = gt;
gtTest.segId = gt.segId(testIds,:);
gtTest.label = gt.label(testIds,:);
gtTest.featMat = gt.featMat(testIds,:);

% train with increasing training data sizes
trainFrac = [0.2, 0.4, 0.6, 0.8, 1];
trainSizes = floor(trainFrac.*length(trainIds));
rng(1); % for reproducibility
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

    [precRec, fig] = TypeEM.Classifier.evaluate(param, curClassifier, gtTest);
    title([methodUsed ' with trainSize:' num2str(curTrainSize)])
    saveas(gcf,fullfile(param.saveFolder,'typeEM',['precrec_test_' methodUsed '_tsize_' num2str(curTrainSize) '.png']))
    close all

    classifiers{end + 1} = curClassifier;
    results(end + 1, :) = {precRec, fig};
end
%{
%% Building output
Util.log('Building output');
classifier = classifiers{end}; %choose trained on all data
Util.save(['/u/sahilloo/H2_3_v2_U1_SubI/type-em/typeClassifier/' datestr(clock,30) '.mat'],classifier,gt,info);
%}
%{
%% Inspect FPs
param.experimentName = 'H2_3_v2_U1_SubI_mr2e_wsmrnet';
gtTestModified = repmat(gtTest,2,1); % because of inverted features
outDir = fullfile(param.saveFolder,'connectEM');

scoreThr = -15.0287; % 75% recall, 93% precision
skel = H2_3_v2_SubI.ConnectEM.inspectFPs(param, gtTestModified, curScores, scoreThr);
skel.write(fullfile(outDir, sprintf('fps-%03f.nml',scoreThr)));
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

