rootDir = '/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet/';
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;
load(fullfile(param.saveFolder,'connectEM','features_all.mat'))
info = Util.runInfo();
Util.showRunInfo(info);

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

classifiers = cell(0);
results = cell(0, 2);

for curTrainSize = trainSizes
       curTrainIds = curRandIds(1:curTrainSize);

        curTrainFeats = [ ...
           [rawFeatsTrain(curTrainIds, :); ...
            fm.invertDirection(rawFeatsTrain(curTrainIds, :))], ...
           [classFeatsTrain(curTrainIds, :); ...
            fm.invertDirection(classFeatsTrain(curTrainIds, :))]];
        curTrainLabels = repmat(gtTrain.label(curTrainIds), 2, 1);

        display(['train '])
        curTrainSize
        size(curTrainFeats)
        size(curTrainLabels)

        
        
        Util.log(' create a neural network')
        net = feedforwardnet(6);
        % set early stopping parameters
        %net.divideParam.trainRatio = 0.85; % training set [%]
        %net.divideParam.valRatio = 0.15; % validation set [%]
        % net.divideParam.testRatio = 0.15; % test set [%]
        % number of hidden layer neurons
        for i=1:net.numLayers-1
            % hidden layer transfer function
            net.layers{i}.transferFcn = 'logsig';
        end
        % network training
        net.trainFcn = 'trainb';
        net.trainParam.epochs = 1000;
        net.performFcn = 'mse';
        
        [net,tr,Y,E]= train(net,curTrainFeats',curTrainLabels');

        curScores = net(curTestFeats');
        
        classifiers{end + 1} = net;
        results(end + 1, :) = {curTestLabels, curScores}; %#ok
        clear net
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
%savefig(fullfile(param.saveFolder,'connectEM',['precrec_test_FFN.fig']))
saveas(gcf,fullfile(param.saveFolder,'connectEM',['precrec_test_FFN.png']))
%Util.save(fullfile(param.saveFolder,'connectEM','evaluateWithFFn.mat'),classifiers,info)


