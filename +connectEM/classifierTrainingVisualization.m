%% Load data from classifierTraining.m
load /gaba/u/mberning/results/edgeClassifier/20170729T131244_workspace.mat;
clear a ans b gt idxTest idxTrain;

%% Calculate predicition values for current classifier
[~, train.scoresClass] = classifierNewFeatures.predict( ...
    cat(1,cat(2, train.rawFeatures, train.classFeatures),cat(2, train.rawFeaturesInv, train.classFeaturesInv)));

[~, test.scoresClass] = classifierNewFeatures.predict( ...
    cat(1,cat(2, test.rawFeatures, test.classFeatures),cat(2, test.rawFeaturesInv, test.classFeaturesInv)));
clear classifierNewFeatures classifier;

%% ... and for earlier classifier with bug discovered by Alessandro
old = load('/gaba/u/mberning/results/edgeClassifier/20170322T153247.mat');

% Note that earlier P/R curves would also have predicted on buggy features,
% NOT reproduced here, e.g. this is actual performance on whole dataset
% using predictDataset.m 
[~, train.scoresClassOld] = old.classifier.predict( ...
    cat(1,cat(2, train.rawFeatures, train.classFeatures),cat(2, train.rawFeaturesInv, train.classFeaturesInv)));

[~, test.scoresClassOld] = old.classifier.predict( ...
    cat(1,cat(2, test.rawFeatures, test.classFeatures),cat(2, test.rawFeaturesInv, test.classFeaturesInv)));
clear old;

%% Transfer to proabilities
sigmoid = @(x)1./(1+exp(-1.*x));

train.probClass = sigmoid(train.scoresClass(:,1));
test.probClass = sigmoid(test.scoresClass(:,1));
train.probClassOld = sigmoid(train.scoresClassOld(:,1));
test.probClassOld = sigmoid(test.scoresClassOld(:,1));

save('/home/mberning/localStorage/classifierTrainingVisualization.mat', '-v7.3');

%% Add (smaller) segment and border size for each edge

borderMeta = load([p.saveFolder 'globalBorder.mat']);
graph = load([p.saveFolder 'graphNewNew.mat']);

% Use only single border edges from test set
[~, idx] = unique(test.edges, 'rows');
uniqueIdx = unique(idx);
counts = histc(idx, uniqueIdx);
idx = uniqueIdx(counts == 1);
fieldNames = fieldnames(test);
for j=1:length(fieldNames)
    test.(fieldNames{j}) = test.(fieldNames{j})(idx,:);
end
clear uniqueIdx idx counts fieldNames;

% Add border size for each test border
[idx, loc] = ismember(graph.edges, test.edges, 'rows');
test.borderSize(loc(idx)) = borderMeta.borderSize(graph.borderIdx(idx));
test.borderSize = test.borderSize';

% Add smaller segment size for each test border
test.smallerSegmentSize = min(segMeta.voxelCount(test.edges),[],2);

%% Load agglo segment predictions
segmentClass = load([p.saveFolder 'segmentAggloPredictions.mat']);
% Extract axon probability
segmentClass.axonProb = zeros(segMeta.maxSegId, 1);
idx = ~isnan(segmentClass.probsMulti(:,2));
segmentClass.axonProb(segmentClass.segId(idx)) = segmentClass.probsMulti(idx,2);
test.axonProb = arrayfun(@(x)segmentClass.axonProb(x), test.edges);
% Extract dendrite probability
segmentClass.dendriteProb = zeros(segMeta.maxSegId, 1);
idx = ~isnan(segmentClass.probsMulti(:,3));
segmentClass.dendriteProb(segmentClass.segId(idx)) = segmentClass.probsMulti(idx,3);
test.dendriteProb = arrayfun(@(x)segmentClass.dendriteProb(x), test.edges);

%% Plot P/R curve
label = [];
figure;
hold on;
set(gca,'FontSize',14);
set(gca,'LineWidth',2);

idx = test.labels ~= 0;
[Xpr,Ypr,Tpr,AUCpr] = perfcurve(test.labels(idx), ...
    test.probClass(idx), ...
    1, 'xCrit', 'reca', 'yCrit', 'prec', 'TVals', 0:0.001:1);
Ypr(1) = 1;
plot(Xpr,Ypr, '-b', 'LineWidth', 2);
label{1} = ['new, on whole test set' ...
    ' (AUC: ' num2str(AUCpr) ')'];

idx = test.labels ~= 0;
[Xpr,Ypr,Tpr,AUCpr] = perfcurve(test.labels(idx), ...
    test.probClassOld(idx), ...
    1, 'xCrit', 'reca', 'yCrit', 'prec', 'TVals', 0:0.001:1);
Ypr(1) = 1;
plot(Xpr,Ypr, '--b', 'LineWidth', 2);
label{2} = ['old, on whole test set' ...
    ' (AUC: ' num2str(AUCpr) ')'];

% Restricted based on segment and border size
idx = test.smallerSegmentSize > 1000 & test.borderSize > 100 & test.labels ~= 0;
[Xpr,Ypr,Tpr,AUCpr] = perfcurve(test.labels(idx), ...
    test.probClass(idx), ...
    1, 'xCrit', 'reca', 'yCrit', 'prec', 'TVals', 0:0.001:1);
Ypr(1) = 1;
plot(Xpr,Ypr, '-c', 'LineWidth', 2);
label{3} = ['new, restricted to above 100 voxel border and 1000 voxel segment size,' ...
    ' (AUC: ' num2str(AUCpr) ')'];

idx = test.smallerSegmentSize > 1000 & test.borderSize > 100 & test.labels ~= 0;
[Xpr,Ypr,Tpr,AUCpr] = perfcurve(test.labels(idx), ...
    test.probClassOld(idx), ...
    1, 'xCrit', 'reca', 'yCrit', 'prec', 'TVals', 0:0.001:1);
Ypr(1) = 1;
plot(Xpr,Ypr, '--c', 'LineWidth', 2);
label{4} = ['old, restricted to above 100 voxel border and 1000 voxel segment size,' ...
    ' (AUC: ' num2str(AUCpr) ')'];

% Additionally restricted to axon segments
idx = test.smallerSegmentSize > 1000 & test.borderSize > 100 & test.labels ~= 0 & min(test.axonProb, [], 2) > 0.3;
[Xpr,Ypr,Tpr,AUCpr] = perfcurve(test.labels(idx), ...
    test.probClass(idx), ...
    1, 'xCrit', 'reca', 'yCrit', 'prec', 'TVals', 0:0.001:1);
Ypr(1) = 1;
plot(Xpr,Ypr, '-g', 'LineWidth', 2);
label{5} = ['new, and only on axon segments recovered 30% probability,' ...
    ' (AUC: ' num2str(AUCpr) ')'];

idx = test.smallerSegmentSize > 1000 & test.borderSize > 100 & test.labels ~= 0 & min(test.axonProb, [], 2) > 0.3;
[Xpr,Ypr,Tpr,AUCpr] = perfcurve(test.labels(idx), ...
    test.probClassOld(idx), ...
    1, 'xCrit', 'reca', 'yCrit', 'prec', 'TVals', 0:0.001:1);
Ypr(1) = 1;
plot(Xpr,Ypr, '--g', 'LineWidth', 2);
label{6} = ['old, and only on axon segments recovered 30% probability,' ...
    ' (AUC: ' num2str(AUCpr) ')'];

% Alternatively restricted to dendrite segments
idx = test.smallerSegmentSize > 1000 & test.borderSize > 100 & test.labels ~= 0 & min(test.dendriteProb, [], 2) > 0.3;
[Xpr,Ypr,Tpr,AUCpr] = perfcurve(test.labels(idx), ...
    test.probClass(idx), ...
    1, 'xCrit', 'reca', 'yCrit', 'prec', 'TVals', 0:0.001:1);
Ypr(1) = 1;
plot(Xpr,Ypr, '-r', 'LineWidth', 2);
label{7} = ['new, and only on dendrite segments recovered 30% probability, ' ...
    ' (AUC: ' num2str(AUCpr) ')'];

idx = test.smallerSegmentSize > 1000 & test.borderSize > 100 & test.labels ~= 0 & min(test.dendriteProb, [], 2) > 0.3;
[Xpr,Ypr,Tpr,AUCpr] = perfcurve(test.labels(idx), ...
    test.probClassOld(idx), ...
    1, 'xCrit', 'reca', 'yCrit', 'prec', 'TVals', 0:0.001:1);
Ypr(1) = 1;
plot(Xpr,Ypr, '--r', 'LineWidth', 2);
label{8} = ['old, and only on dendrite segments recovered 30% probability, ' ...
    ' (AUC: ' num2str(AUCpr) ')'];

xlim([0 1]);
ylim([0 1]);
axis square;
legend(label, 'Location', 'Best');
