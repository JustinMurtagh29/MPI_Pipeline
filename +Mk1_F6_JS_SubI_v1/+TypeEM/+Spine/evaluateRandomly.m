% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
% Modified by
%   Sahil Loomba <sahil.loomba@brain.mpg.de>

% Randomly shuffled all boxes into training and test set 80-20 %

clear;
timeStamp = datestr(now,'yyyymmddTHHMMSS');
methodUsed = 'LogitBoost'; %'AdaBoostM1'; % 'LogitBoost';
numTrees = 1500;
addpath(genpath('/gaba/u/sahilloo/repos/amotta/matlab/'))

rootDir = '/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet/';
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;
param.experimentName = 'Mk1_F6_JS_SubI_v1_mrnet_wsmrnet';
import Mk1_F6_JS_SubI_v1.TypeEM.*

info = Util.runInfo();

segmentMeta = load([param.saveFolder 'segmentMeta.mat'], 'voxelCount', 'point');
vxThr = 50;

% load training data
featureSetName = 'segmentAgglomerate'; %'segmentAgglomerate'; % 'segment'
Util.log(sprintf('Evaluating for %s features',featureSetName))

% load spinehead training data
nmlDir = fullfile(param.saveFolder, ...
     'tracings', 'box-seeded','spine-head-ground-truth');
nmlFiles = fullfile(nmlDir, ...
     {'spine-head-ground-truth-1.nml','spine-head-ground-truth-2.nml', 'spine-head-ground-truth-3.nml', ...
     'spine-head-ground-truth-4.nml','spine-head-ground-truth-5.nml','spine-head-ground-truth-6.nml',...
     'spine-head-ground-truth-7.nml','spine-head-ground-truth-8.nml','spine-head-ground-truth-9.nml',...
     'spine-head-ground-truth-10.nml','spine-head-ground-truth-11.nml',...
     'spine-head-ground-truth-13.nml', ...
     'spine-head-ground-truth-14.nml','spine-head-ground-truth-15.nml','spine-head-ground-truth-16.nml',...
     'spine-head-ground-truth-18.nml','spine-head-ground-truth-19.nml',...
     'spine-head-ground-truth-20.nml', 'spine-head-ground-truth-21.nml','spine-head-ground-truth-23.nml',...
    'spine-head-ground-truth-24.nml','spine-head-ground-truth-25.nml'});

rng(0);
idxTrain = 1:numel(nmlFiles); %[1,2,3,4,5,6,7,8,9,10,11];
% load train set
curNodes = table();
gt = struct;
gt.class = categorical({'spinehead'});
gt.segId = [];
gt.label = [];
curMask = gt.class == 'spinehead';
for i = idxTrain
    curNml = slurpNml(nmlFiles{i});
    curNodes = NML.buildNodeTable(curNml);
    
    curNodes.coord = curNodes.coord + 1;
    curNodes.segId = Seg.Global.getSegIds(param, curNodes.coord);
    assert(all(curNodes.segId));
    
    curBox = curNml.parameters.userBoundingBox;
    curBox = { ...
        curBox.topLeftX, curBox.topLeftY, curBox.topLeftZ, ...
        curBox.width, curBox.height, curBox.depth};
    curBox = cellfun(@str2double, curBox);
    
    curBox = Util.convertWebknossosToMatlabBbox(curBox);
    curseg = loadSegDataGlobal(param.seg, curBox);
    
    posSegIds = reshape(unique(curNodes.segId), [], 1);
    negSegIds = reshape(setdiff(curseg, [0; posSegIds]), [], 1);

    % throw away segments that are too small
    vxCount = segmentMeta.voxelCount(negSegIds); 
    toDel = vxCount < vxThr;
    negSegIds(toDel) = [];

    gt.segId = cat(1,gt.segId, double(posSegIds), double(negSegIds));
    labelTemp = zeros(numel(posSegIds)+numel(negSegIds), numel(gt.class));
    labelTemp(1:numel(posSegIds), curMask) = +1;
    labelTemp((numel(posSegIds) + 1):end, curMask) = -1;
    
    gt.label = cat(1,gt.label, labelTemp);

    clear curNml
end
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

classifiers = cell(0);
results = cell(0, 2);

for curTrainSize = trainSizes
    curTrainIds = trainIds(1:curTrainSize);

    gtTrain = gt;
    gtTrain.segId = gt.segId(curTrainIds,:);
    gtTrain.label = gt.label(curTrainIds,:);
    gtTrain.featMat = gt.featMat(curTrainIds,:);

    curClassifier.classes = gt.class; % classifier for spinehead class only
    curClassifier.classifiers = arrayfun( ...
            @(c) buildForClass(gtTrain, c, methodUsed, numTrees), ...
            gtTrain.class, 'UniformOutput', false);
    curClassifier.featureSetName = featureSetName;

    % apply classifier to test data
    [precRec, fig, curGtTest] = TypeEM.Classifier.evaluate(param, curClassifier, gtTest);
    title([methodUsed ' with trainSize:' num2str(curTrainSize) '_trees:' num2str(numTrees)])
    saveas(gcf,fullfile(param.saveFolder,'typeEM','spine',featureSetName,[timeStamp '_precrec_box_' methodUsed '_tsize_' num2str(curTrainSize) '_trees_' num2str(numTrees) '_Random.png']))
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

% Look at false positives
curCount = 40;
className = 'spinehead';
skels = Debug.inspectFPs(param, curCount, className, curGtTest);
skels.write(fullfile(param.saveFolder,'typeEM','spine', featureSetName,sprintf('%s_fps-%s-Random.nml',timeStamp, className)));

% Look at true positives
curCount = 40;
className = 'spinehead';
skels = Debug.inspectTPs(param, curCount, className, curGtTest);
skels.write(fullfile(param.saveFolder,'typeEM','spine', featureSetName,sprintf('%s_tps-%s-Random.nml',timeStamp, className)));

% label statistics
Spine.Debug.labelDistributions

%{

%% Building output
Util.log('Building output');
classifier = classifiers{end}; %choose trained on all data
Util.save(['/u/sahilloo/Mk1_F6_JS_SubI_v1/type-em/spineClassifier/' datestr(clock,30) '.mat'],classifier,gt,info);
%}

function classifier = buildForClass(gt, class, methodUsed, numTrees)
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
        methodUsed, numTrees, 'tree', ...
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

