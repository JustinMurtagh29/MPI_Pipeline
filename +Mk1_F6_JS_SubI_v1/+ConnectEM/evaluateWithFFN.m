rootDir = '/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet/';
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;
load(fullfile(param.saveFolder,'connectEM','features_all.mat'))
info = Util.runInfo();
Util.showRunInfo(info);

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
gtTestModified = repmat(gtTest,2,1); % because of inverted features

% compile train data
curTrainFeats = [ ...
           [rawFeatsTrain; ...
            fm.invertDirection(rawFeatsTrain)], ...
           [classFeatsTrain; ...
            fm.invertDirection(classFeatsTrain)]];
curTrainLabels = repmat(gtTrain.label, 2, 1);

Util.log(' create a neural network')
net = feedforwardnet(6);
% set early stopping parameters
net.divideParam.trainRatio = 0.85; % training set [%]
net.divideParam.valRatio = 0.15; % validation set [%]
% % net.divideParam.testRatio = 0.15; % test set [%]

%% number of hidden layer neurons
for i=1:net.numLayers-1
    % hidden layer transfer function
    net.layers{i}.transferFcn = 'logsig';
end

%net = configure(net,inputs',outputs');

%% network training
net.trainFcn = 'trainb';
net.performFcn = 'mse';

[net,tr,Y,E]= train(net,curTrainFeats',curTrainLabels');

%% network response after training
curTestScores = net(curTestFeats');

%% Plot aggregate results
curFig = figure();
curFig.Color = 'white';

curAx = axes(curFig);
curAx.TickDir = 'out';
hold(curAx, 'on');
    
[curPrec, curRec, threshVec] = ...
    TypeEM.Classifier.buildPrecRec(curTestScores, curTestLabels);
curLine = ...
    TypeEM.Classifier.plotPrecisionRecall(curAx, curPrec, curRec);

curLines = flip(curAx.Children);
axis(curAx, 'square');

title(curAx, ...
    {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
savefig(fullfile(param.saveFolder,'connectEM',['precrec_test_FFN.fig']))
saveas(gcf,fullfile(param.saveFolder,'connectEM',['precrec_test_FFN.png']))


