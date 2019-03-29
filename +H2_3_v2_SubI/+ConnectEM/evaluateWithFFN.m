clear all
addpath(genpath('/gaba/u/sahilloo/repos/amotta/'));

rootDir = '/tmpscratch/sahilloo/data/H2_3_v2_U1_SubI/pipelineRun_mr2e_wsmrnet/';
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;
load(fullfile(param.saveFolder,'connectEM','features_all.mat'))
info = Util.runInfo();
Util.showRunInfo(info);

layers = [6]
trainFunc = 'traingdm'
epochs = 2000

%% data
Util.log('divide into training and test set')
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

%{% compile train data
curTrainFeats = [ ...
           [rawFeatsTrain; ...
            fm.invertDirection(rawFeatsTrain)], ...
           [classFeatsTrain; ...
            fm.invertDirection(classFeatsTrain)]];
curTrainLabels = repmat(gtTrain.label, 2, 1);
%}

% train with increasing training data sizes
trainFrac = [0.2, 0.4, 0.6, 0.8, 1];
trainSizes = floor(trainFrac.*height(gtTrain));
rng(1); % for reproducibility
curRandIds = randperm(height(gtTrain));

classifiers = cell(0,2);
results = cell(0, 2);
%% Plot aggregate results
curFig = figure();
curFig.Color = 'white';

curAx = axes(curFig);
curAx.TickDir = 'out';
hold(curAx, 'on');

for curTrainSize = trainSizes
       curTrainIds = curRandIds(1:curTrainSize);

        curTrainFeats = [ ...
           [rawFeatsTrain(curTrainIds, :); ...
            fm.invertDirection(rawFeatsTrain(curTrainIds, :))], ...
           [classFeatsTrain(curTrainIds, :); ...
            fm.invertDirection(classFeatsTrain(curTrainIds, :))]];
        curTrainLabels = repmat(gtTrain.label(curTrainIds), 2, 1);
        % train and predict with FFN
        [curScores, net, tr] = trainFFN(curTrainFeats, curTrainLabels, curTestFeats,...
                             layers, trainFunc, epochs);

        [curPrec, curRec, threshVec] = ...
            TypeEM.Classifier.buildPrecRec(curScores, curTestLabels);
        curLine = ...
            TypeEM.Classifier.plotPrecisionRecall(curAx, curPrec, curRec);  

        classifiers{end + 1} = {net, tr};        
        results(end + 1, :) = {curTestLabels, curScores}; %#ok
end

curLines = flip(curAx.Children);
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
%savefig(fullfile(param.saveFolder,'connectEM',['precrec_FFN.fig']))
saveas(gcf,fullfile(param.saveFolder,'connectEM',['precrec_' trainFunc '_'...
                     num2str(layers) '_' num2str(epochs) '.png']))
Util.save(fullfile(param.saveFolder,'connectEM',['data_' trainFunc '_'...
                     num2str(layers) '_' num2str(epochs) '.mat']),classifiers,info)

function [curScores,net, tr] = trainFFN(curTrainFeats, curTrainLabels, curTestFeats,...
                                         layers, trainFunc, epochs)
    net = feedforwardnet(layers,trainFunc);
    % set early stopping parameters
    net.divideParam.trainRatio = 1; % training set [%]
    net.divideParam.valRatio = 0; % validation set [%]
    net.divideParam.testRatio = 0; % test set [%]
    % number of hidden layer neurons
    for i=1:net.numLayers-1
        net.layers{i}.transferFcn = 'logsig';
    end
    % network training
    %net.trainFcn = 'traingdx'; % 'trainb';
    net.trainParam.epochs = epochs;
    net.performFcn = 'mse';
    
    net = configure(net,curTrainFeats', curTrainLabels');
    %net = init(net);
    [net1,tr,Y,E]= train(net,curTrainFeats',curTrainLabels');
    curScores = net1(curTestFeats');
end
