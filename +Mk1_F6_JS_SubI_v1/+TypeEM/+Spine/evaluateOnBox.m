% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
% Modified by
%   Sahil Loomba <sahil.loomba@brain.mpg.de>

% training set
% Positive labels: 11 boxes dense-spine-head labels
% Negative lables: 13 typeEM boxes with Axon/Dend/Glia labels
% test set:
% Positive labels: All 20 SH nodes in one box
% Negative labels: All the rest segments in the box

clear;
methodUsed = 'LogitBoost'; %'AdaBoostM1'; % 'LogitBoost';
addpath(genpath('/gaba/u/sahilloo/repos/amotta/matlab/'))

timeStamp = datestr(now,'yyyymmddTHHMMSS');
rootDir = '/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet/';
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;
param.experimentName = 'Mk1_F6_JS_SubI_v1_mrnet_wsmrnet';
import Mk1_F6_JS_SubI_v1.TypeEM.*

info = Util.runInfo();
segmentMeta = load([param.saveFolder 'segmentMeta.mat'], 'voxelCount', 'point');
vxThr = 100;

% load training data
featureSetName = 'segmentAgglomerate';
Util.log(sprintf('Evaluating for %s features',featureSetName))
nmlDir = fullfile(param.saveFolder, ...
     'tracings', 'typeEM');
nmlFiles = fullfile(nmlDir, 'proofread', ...
     {'box-1.nml','box-2.nml', 'box-3.nml', ...
     'box-4.nml','box-5.nml','box-6.nml',...
     'box-7.nml','box-8.nml','box-9.nml',...
     'box-10.nml','box-11.nml','box-12.nml',...
     'box-13.nml'});

rng(0);
idxTrain = [1,2,3,4,5,6,7,8,9,10,11,12,13];
% load features
gtType = TypeEM.GroundTruth.loadSet( ...
        param, featureSetName, nmlFiles(idxTrain));
gtType = TypeEM.GroundTruth.loadFeatures(param, featureSetName, gtType);

%{
% remove dendrite segments from training data of typeEM boxes
idxDel = gtType.label(:,3) == -1;
gtType.segId(idxDel,:) = [];
gtType.label(idxDel,:) = [];
gtType.featMat(idxDel,:) = [];
%}
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
gtFiles = randperm(numel(nmlFiles));
idxTrain = gtFiles(1:end-1);
idxTest = gtFiles(end);

% load train set
curNodes = table();
for i = idxTrain
    curNml = slurpNml(nmlFiles{i});
    curNodes = cat(1, curNodes, NML.buildNodeTable(curNml));
    clear curNml
end
curNodes.coord = curNodes.coord + 1;
curNodes.segId = Seg.Global.getSegIds(param, curNodes.coord);
assert(all(curNodes.segId));

gtSH = struct;
gtSH.class = categorical({'glia', 'axon', 'dendrite', 'spinehead'});
gtSH.segId = reshape(unique(curNodes.segId), [], 1);
gtSH.label = [-1*ones(numel(gtSH.segId),3), ones(numel(gtSH.segId),1)]; %[-1,-1,-1,+1]
gtSH = TypeEM.GroundTruth.loadFeatures(param, featureSetName, gtSH);

% combine type gt and spinehead gt
gt = struct;
gt.segId = cat(1,gtType.segId, gtSH.segId);
gt.class = categorical({'glia', 'axon', 'dendrite', 'spinehead'});

gt.label = zeros(numel(gt.segId), numel(gt.class));
gt.label = cat(1, [gtType.label, -1*ones(numel(gtType.segId),1)],... % [x,x,x,-1]
                  gtSH.label); % [-1,-1,-1,+1] 

gt.featMat = cat(1, gtType.featMat, gtSH.featMat); 
gt.featNames = cat(2, gtType.featNames, gtSH.featNames);

% throw away segments that are too small
vxCount = segmentMeta.voxelCount(gt.segId);
toDel = vxCount < vxThr;
gt.segId(toDel,:) = [];
gt.label(toDel,:) = [];
gt.featMat(toDel,:) = [];

% Reduce to spinehead class only
curMask = gt.class == 'spinehead';
gt.class = gt.class(curMask);
gt.label = gt.label(:,curMask);

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
% throw away segments that are too small
vxCount = segmentMeta.voxelCount(negSegIds);
toDel = vxCount < vxThr;
negSegIds(toDel) = [];

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

    curClassifier.classes = gt.class; % classifier for spinehead class only
    curClassifier.classifiers = arrayfun( ...
            @(c) buildForClass(gtTrain, c, methodUsed), ...
            gtTrain.class, 'UniformOutput', false);
    curClassifier.featureSetName = featureSetName;

    % apply classifier to test data
    [precRec, fig, curGtTest] = TypeEM.Classifier.evaluate(param, curClassifier, gtTest);
    title([methodUsed ' with trainSize:' num2str(curTrainSize)])
    saveas(gcf,fullfile(param.saveFolder,'typeEM','spine',featureSetName,[timeStamp '_precrec_box_' methodUsed '_tsize_' num2str(curTrainSize) '_typeEMBoxWDend.png']))
    close all
%{
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
%}
end 

%{
% Look at false positives
curCount = 40;
className = 'spinehead';
skels = Debug.inspectFPs(param, curCount, className, curGtTest);
skels.write(fullfile(param.saveFolder,'typeEM','spine', featureSetName,sprintf('%s-fps-%s_typeEMBox.nml',timeStamp, className)));

% Look at true positives
curCount = 40;
className = 'spinehead';
skels = Debug.inspectTPs(param, curCount, className, curGtTest);
skels.write(fullfile(param.saveFolder,'typeEM','spine', featureSetName,sprintf('%s-tps-%s_typeEMBox.nml',timeStamp, className)));

% label statistics
Spine.Debug.labelDistributions
%}
%{

%% Building output
Util.log('Building output');
classifier = classifiers{end}; %choose trained on all data
Util.save(['/u/sahilloo/Mk1_F6_JS_SubI_v1/type-em/spineClassifier/' datestr(clock,30) '.mat'],classifier,gt,info);
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

