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
param.experimentName = 'Mk1_F6_JS_SubI_v1_mr2e_wsmrnet';

info = Util.runInfo();

% load training data
featureSetName = 'segmentAgglomerate';
nmlDir = fullfile(param.saveFolder, ...
     'tracings', 'typeEM');
nmlFiles = fullfile(nmlDir, 'proofread', ...
     {'box-1.nml','box-2.nml', 'box-3.nml', ...
     'box-4.nml','box-5.nml','box-6.nml',...
     'box-7.nml','box-8.nml','box-9.nml'});

rng(0);
idxTrain = [1,2,3,4,5,6,7,8];
idxTest = 9;
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

    [precRec, fig, curGtTest] = TypeEM.Classifier.evaluate(param, curClassifier, gtTest);
    title([methodUsed ' with trainSize:' num2str(curTrainSize)])
    saveas(gcf,fullfile(param.saveFolder,'typeEM',['precrec_box_' methodUsed '_tsize_' num2str(curTrainSize) '.png']))
    close all

    classifiers{end + 1} = curClassifier;
    results(end + 1, :) = {precRec, fig};
end 

% Look at false positives
className = 'glia';
clear cur*;
curCount = 50;
curDigits = ceil(log10(1 + curCount));
curPoints = Seg.Global.getSegToPointMap(param);

curMask = curGtTest.class == className;

fp = table;
fp.segId = curGtTest.segId;
fp.label = curGtTest.label(:, curMask);
fp.score = curGtTest.scores(:, curMask);

fp = fp(fp.label < 0, :);
fp = sortrows(fp, 'score', 'descend');

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = Skeleton.setDescriptionFromRunInfo(skel, info);

curSegIds = fp.segId(1:curCount);
curScores = fp.score(1:curCount);

curNodes = curPoints(curSegIds, :);
curNames = arrayfun( ...
    @(idx, segId, score) sprintf( ...
        '%0*d. Segment %d. Score %.3f', ...
        curDigits, idx, segId, score), ...
    (1:curCount)', curSegIds, curScores, ...
    'UniformOutput', false);

skel = skel.addNodesAsTrees(curNodes, curNames);
skel.write(fullfile(param.saveFolder,'typeEM', sprintf('fps-%s-%03f.nml',className,scoreThr)));


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

