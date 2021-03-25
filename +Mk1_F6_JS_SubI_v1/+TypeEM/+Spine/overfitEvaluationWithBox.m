% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
% Modified by
%   Sahil Loomba <sahil.loomba@brain.mpg.de>

% training set
% Positive labels: 10 boxes dense-spine-head labels
% Negative lables: all remaining segments in the boxes
% test set:
% Positive labels: All 20 SH nodes in one box
% Negative labels: All the rest segments in the box

% oversample + labels in gtTrain
% undersample - labels in gtTrain

clear;
timeStamp = [datestr(now,'yyyymmddTHHMMSS')];
methodUsed = 'LogitBoost'; %'AdaBoostM1'; % 'LogitBoost';
numTrees = 1500;
addNoise = true;
noiseDev = 0.01;
evalForAgglos = true;
vxThrTest = true;
gtVersion = 'v4';
featureSetName = 'segmentAgglomerate'; %'segmentAgglomerate'; % 'segment'
factorPos = 0.65; % 100x times oversample
factorNeg = 0.8; % x times undersample
featsImp = false;

rootDir = '/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet';
addpath(genpath('/gaba/u/sahilloo/repos/amotta/matlab/'))

param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;
param.experimentName = 'Mk1_F6_JS_SubI_v1_mrnet_wsmrnet';
maxSegId = Seg.Global.getMaxSegId(param);
import Mk1_F6_JS_SubI_v1.TypeEM.*

info = Util.runInfo();

segmentMeta = load([param.saveFolder 'segmentMeta.mat'], 'voxelCount', 'point');
vxThr = 100;

if featsImp
    load(fullfile('/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet/typeEM/spine/segmentAgglomerate/importantFeatures.mat'));
end

Util.log(sprintf('Evaluating for %s features',featureSetName))

%for idxTest = gtFiles
    % load spinehead training data
    nmlDir = fullfile(param.saveFolder, ...
         'tracings', 'box-seeded','spine-head-ground-truth', gtVersion);
    nmlFiles = arrayfun(@(x) fullfile(nmlDir, x.name), dir(fullfile(nmlDir, '*.nml')), 'uni',0);
    
    rng(0); % default 0
    gtFiles = 1:numel(nmlFiles);
    idxTest = 39; %gtFiles(end);
    idxTrain = gtFiles ;%setdiff(gtFiles, idxTest); % all except test nml
 
    if addNoise
        experimentName = sprintf('overfit_box_%s_%d_over_%.2f_under_%.2f_cost_100_addNoise_%.3f_testset_%d_GT_%s_features_%s_evalForAgglo_%d', ...
                methodUsed, numTrees, factorPos, factorNeg, noiseDev, idxTest, gtVersion, featureSetName, evalForAgglos);
    else
        experimentName = sprintf('overfit_box_%s_%d_over_%.2f_under_%.2f_cost_100_noNoise_testset_%d_GT_%s_features_%s_evalForAgglo_%d', ...
            methodUsed, numTrees, factorPos, factorNeg, idxTest, gtVersion, featureSetName, evalForAgglos);
    end
    
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
    
        %oversample pos segIds
        posSegIds = repmat(posSegIds,int32(100*factorPos),1);
        negSegIds = reshape(setdiff(curseg, [0; posSegIds]), [], 1);
        
        % undersample negative labels
        negSegIds = negSegIds(randperm(length(negSegIds), ceil(factorNeg*length(negSegIds)))); % choose x% of ids
    
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
    
    % load test set
    curNml = slurpNml(nmlFiles{idxTest});
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
    
    if vxThrTest
        % throw away segments that are too small
        vxCount = segmentMeta.voxelCount(negSegIds);
        toDel = vxCount < vxThr;
        negSegIds(toDel) = [];
    end

    gtTest = struct;
    gtTest.class = categorical({'spinehead'});
    gtTest.segId = cat(1, double(posSegIds), double(negSegIds));
    gtTest.label = zeros(numel(gtTest.segId), numel(gtTest.class));
    
    curMask = gtTest.class == 'spinehead';
    gtTest.label(1:numel(posSegIds), curMask) = +1;
    gtTest.label((numel(posSegIds) + 1):end, curMask) = -1;
    
    curRandIds = randperm(size(gt.label,1));
    trainIds = curRandIds;
    
    % train with increasing training data sizes
    trainFrac = 1; %[0.2, 0.4, 0.6, 0.8, 1];
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
    
        if addNoise
            gtTrain.featMat = augmentFeatures(gtTrain.featMat, noiseDev);
        end
    
        if featsImp
            gtTrain.featMat = gtTrain.featMat(:,idxImpFeats);
            gtTrain.featNames = gtTrain.featNames(:,idxImpFeats);
        end
        curClassifier.classes = gt.class; % classifier for spinehead class only
        curClassifier.classifiers  = arrayfun( ...
                @(c) buildForClass(gtTrain, c, methodUsed, numTrees), ...
                gtTrain.class, 'UniformOutput', false);
        curClassifier.featureSetName = featureSetName;
    
        % apply classifier to test data
        if evalForAgglos
            % add agglos to gtTest to update to max of neighbour scores       
            [~,idxSort] = sort(gtTest.segId);
            gtTest.segId = gtTest.segId(idxSort);
            gtTest.label = gtTest.label(idxSort);
    
            [~,idxSort] = sort(curNodes.segId);
            [~, newTreeIds] = histc(curNodes.treeId(idxSort), unique(curNodes.treeId)); % order them from 1 onwards, sorted by ascendign segIds
            
            shSegIds = gtTest.segId(gtTest.label==1);
            nonShSegIds = gtTest.segId(gtTest.label==-1);
            agglos = zeros(maxSegId,1);

            % curNodes.segIds have to be sorted as well to match shSegIds
            agglos(shSegIds) = newTreeIds;
            agglos(nonShSegIds) = reshape(max(newTreeIds) + (1:numel(nonShSegIds)), 1,'');
            agglos(agglos==0) = [];
            gtTest.agglos = agglos;
            clear agglos;
        end

        if featsImp
            gtTest.featsImp = idxImpFeats;
        end

        [precRec, fig, curGtTest] = TypeEM.Classifier.evaluate(param, curClassifier, gtTest);

        title(sprintf('%s \n trainSize: %d numTrees: %d', experimentName, curTrainSize, numTrees), 'Interpreter', 'none')
        saveas(gcf,fullfile(param.saveFolder,'typeEM','spine',featureSetName,...
                [timeStamp '_precrec_' experimentName '.png']))
        %{
        saveas(gcf,fullfile(param.saveFolder,'typeEM','spine',featureSetName,...
                [timeStamp '_precrec_' experimentName '.fig']))
        %}
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
    
    sprintf('GT train: + labels: %d, - labels: %d', sum(gtTrain.label==1), sum(gtTrain.label==-1))
    sprintf('GT test: + labels: %d, - labels: %d', sum(gtTest.label==1), sum(gtTest.label==-1))
    
    % Look at false positives
    curCount = 100;
    className = 'spinehead';
    skels = Debug.inspectFPs(param, curCount, className, curGtTest);
    skels.write(fullfile(param.saveFolder,'typeEM','spine', featureSetName,sprintf('%s_true-labels-%s-%s.nml',timeStamp, className, experimentName)));
    
    % Look at true positives
    curCount = 100;
    className = 'spinehead';
    skels = Debug.inspectTPs(param, curCount, className, curGtTest);
    skels.write(fullfile(param.saveFolder,'typeEM','spine', featureSetName,sprintf('%s_false-labels-%s-%s.nml',timeStamp, className, experimentName)));
    
%end

%{
% label statistics
Spine.Debug.labelDistributions

%% Building output
Util.log('Building output');
classifier = classifiers{end}; %choose trained on all data
Util.save(['/u/sahilloo/Mk1_F6_JS_SubI_v1/type-em/spineClassifier/' datestr(clock,30) '.mat'],classifier,gt, curGtTest, precRec, info);

% extract the score threshold
scoreVec = curGtTest.scores;
[threshVec, threshIndVec] = ...
        sort(scoreVec, 'descend');

probsVec = curGtTest.probs;
[probVec, probsIndVec] = ...
        sort(probsVec, 'descend');

prec = 1; rec = 0.9;
idxThr = find(precRec.spinehead(:,1)==prec & precRec.spinehead(:,2)==rec);
sprintf('Score thr is %f at pre %f/ rec %f',threshVec(idxThr),prec, rec)
sprintf('Prob thr is %f at pre %f/ rec %f',probVec(idxThr),prec, rec)
%}
function classifier = buildForClass(gt, class, methodUsed, numTrees)
    % classifier = buildForClass(data, className)
    %   Builds a one-versus-all (OVA) classifier for a class

    classIdx = find(gt.class == class);
    mask = gt.label(:, classIdx) ~= 0;

    % extract training data
    trainFeat = gt.featMat(mask, :);
    trainClass = gt.label(mask, classIdx) > 0;

    % weights for + labels
    weights = ones(length(trainClass),1);
    weights(trainClass==1) = 1;

    % cost from confusion matrix
    cost.ClassNames = [true, false];
    cost.ClassificationCosts = [0,100; 1,0];

    % train classifier
    switch methodUsed
        case 'LogitBoost'
            classifier = fitensemble( ...
                trainFeat, trainClass, ...
                methodUsed, numTrees,'tree', ...
                'LearnRate', 0.1, ...
                'NPrint', 100, ...
                'Prior', 'Empirical', ...
                'Type', 'Classification', ...
                'Cost', cost, ...
                'crossval','off', ...
                'ClassNames', [true, false]);
        case 'RobustBoost'
            classifier = fitcensemble( ...
                trainFeat, trainClass, ...
                'method', methodUsed, ...
                'NPrint', 100, ...
                'Prior', 'Empirical', ...
                'Type', 'Classification', ...
                'Learners', 'Tree', ...
                'RobustErrorGoal',0.01, ...
                'NumLearningCycles',numTrees, ...
                'Cost', cost, ...
                'ClassNames', [true, false]);
        case 'GentleBoost'
            classifier = fitcensemble( ...
                trainFeat, trainClass, ...
                'method', methodUsed, ...
                'NPrint', 100, ...
                'LearnRate', 0.1, ...
                'Prior', 'Empirical', ...
                'Type', 'Classification', ...
                'Learners', 'Tree', ...
                'NumLearningCycles',numTrees, ...
                'Cost', cost, ...
                'ClassNames', [true, false]);
        case 'AdaBoostM1'
            classifier = fitcensemble( ...
                trainFeat, trainClass, ...
                'method', methodUsed, ...
                'NPrint', 100, ...
                'LearnRate', 0.84803, ...
                'Prior', 'Empirical', ...
                'Type', 'Classification', ...
                'Learners', 'Tree', ...
                'NumLearningCycles',numTrees, ...
                'Cost', cost, ...
                'ClassNames', [true, false]);
        case 'Optimize'
             classifier = fitcensemble(trainFeat,trainClass,'OptimizeHyperparameters','auto',...
                'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName','expected-improvement-plus'));
    end
end

function [a, b] = trainPlattForClass(gt, classIdx)
    mask = gt.label(:, classIdx) ~= 0;
    scores = gt.scores(mask, classIdx);
    labels = double(gt.label(mask, classIdx) > 0);

    % optimize Platt parameter
    [a, b] = TypeEM.Classifier.learnPlattParams(scores, labels);
end

function X = augmentFeatures(X, dev)
    %AUGMENTFEATURES Add random noise to features based on the feature mean.
    % INPUT X: [NxM] float
    %           Feature matrix. Rows correspond to instances and columns to
    %           features.
    %       dev: (Optional) float
    %           Fraction of the feature mean that is used as the standard
    %           deviation for random noise for the corresponding feature. i.e.
    %           for each column the following is done
    %           X(:,i) = X(:,i) + randn(length(X(:,i), 1), 1).*mean(X(:,i))*dev
    %           (Default: 0.01)
    % OUTPUT X: [NxM] float
    %           The input features with added random noise.
    % Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>
    
    if ~exist('dev', 'var') || isempty(dev)
        dev = 0.01;
    end
    
    m = mean(X, 1);
    X = bsxfun(@plus, X, bsxfun(@times, randn(size(X), 'like', X), m.*dev));
    
end