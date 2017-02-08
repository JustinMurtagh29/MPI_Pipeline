%% Settings
% Load parameter of pipeline run
load /gaba/u/mberning/results/pipeline/20161120_ROI/allParameter.mat;

% Extract information of dense cube annotations
newGTfolder = '+connectEM/trainingData/';
files = dir([newGTfolder '*.nml']);
temp = cellfun(@(x)strsplit(x, '__'), {files(:).name}, 'uni', 0);
taskIds = cellfun(@(x)x{2}, temp, 'uni', 0);
tracer = cellfun(@(x)x{3}, temp, 'uni', 0);
annotationIds = cellfun(@(x)x{4}(1:end-4), temp, 'uni', 0);

% Additional information needed
% Note that trainFile has to fit bboxSmall
files = {files(:).name};
pT.local(1).bboxSmall = [4133 3963 2253; 4578 4408 2432]';
pT.local(2).bboxSmall = [4438 1320 893; 4883 1765 1072]';
pT.local(3).bboxSmall = [1824 6673 1239; 2269 7118 1418]';
pT.local(1).trainFile = cellfun(@(x)[newGTfolder x], files(1:3), 'uni', 0);
pT.local(2).trainFile = cellfun(@(x)[newGTfolder x], files(4:6), 'uni', 0);
pT.local(3).trainFile = cellfun(@(x)[newGTfolder x], files(7:9), 'uni', 0);

%% Collect training data from segmentation + nmls
% Get a list of labels from each (redundant) training file in each region
gt = connectEM.getContinuityLabelsFromNml(p, pT);

% Consolidate redundant annotations
gt = connectEM.consolidateRedundantAnnotations(gt);

%% Collect features from SynEM feature calculation for all edges in GT
gt = connectEM.getFeaturesForEdges(p, gt);

%% Divide GT into training an test set (80%/20%) in each dense cube

% Make dimensionality of labels and prob consitent
for i=1:length(gt)
    gt(i).labels = gt(i).labels';
    gt(i).prob = gt(i).prob';
end

train = struct();
test = struct();
for i=1:length(gt)
    idxTrain = randperm(numel(gt(i).labels), floor(numel(gt(i).labels).*0.8));
    idxTest = setdiff(1:numel(gt(i).labels), idxTrain);
    fieldNames = fieldnames(gt);
    for j=1:length(fieldNames)
        train(i).(fieldNames{j}) = gt(i).(fieldNames{j})(idxTrain,:);
        test(i).(fieldNames{j}) = gt(i).(fieldNames{j})(idxTest,:);
    end
end

train = Util.concatStructs(1, train(1), train(2), train(3));
test = Util.concatStructs(1, test(1), test(2), test(3));

% plot class distribution of some collected features
% for i=1:40
%     figure('Position', [1 1 1920 999]);
%     histogram(train.oldFeatures(train.labels == 1,i))
%     hold on; histogram(train.oldFeatures(train.labels == -1,i))
% end

%% Train different classifier on all 3 feature training sets
classifierOldFeatures = connectEM.trainClassifier( train.oldFeatures, train.labels );
classifierRawFeatures = connectEM.trainClassifier( train.rawFeatures, train.labels );
classifierNewFeatures = connectEM.trainClassifier( cat(2, train.rawFeatures, ...
    train.classFeatures), train.labels );

%% Calculate predicition values for all classifiers (training & test)
[~, train.scoresOld] = classifierOldFeatures.predict(train.oldFeatures);
[~, train.scoresRaw] = classifierRawFeatures.predict(train.rawFeatures);
[~, train.scoresClass] = classifierNewFeatures.predict(cat(2, train.rawFeatures, ...
    train.classFeatures));

[~, test.scoresOld] = classifierOldFeatures.predict(test.oldFeatures);
[~, test.scoresRaw] = classifierRawFeatures.predict(test.rawFeatures);
[~, test.scoresClass] = classifierNewFeatures.predict(cat(2, test.rawFeatures, ...
    test.classFeatures));

%% Transfer to proabilities
sigmoid = @(x)1./(1+exp(-1.*x));

train.probOld = sigmoid(train.scoresOld(:,1));
train.probRaw = sigmoid(train.scoresRaw(:,1));
train.probClass = sigmoid(train.scoresClass(:,1));

test.probOld = sigmoid(test.scoresOld(:,1));
test.probRaw = sigmoid(test.scoresRaw(:,1));
test.probClass = sigmoid(test.scoresClass(:,1));


%% Visualize predicted probabilities vs. frequencies in test set
binSize = 500;
binLowerLimit = 1:binSize:length(test.prob);
binUpperLimit = [(binSize):binSize:length(test.prob) length(test.prob)];
binCenter = (binLowerLimit + binUpperLimit) ./ 2;
figure;
subplot(2,2,1);
[sortedProb, idx] = sort(test.prob);
sortedLabels = test.labels(idx);
sortedProbBinned = arrayfun(@(x,y)sum(sortedLabels(x:y) == 1)./numel(sortedLabels(x:y)), ...
    binLowerLimit, binUpperLimit);
plot(sortedProb, '-k');
hold on;
plot(binCenter, sortedProbBinned, 'xk');
title('Old classifier, old features: Probability vs. Frequency in test set');
subplot(2,2,2);
[sortedProb, idx] = sort(test.probOld);
sortedLabels = test.labels(idx);
sortedProbBinned = arrayfun(@(x,y)sum(sortedLabels(x:y) == 1)./numel(sortedLabels(x:y)), ...
    binLowerLimit, binUpperLimit);
plot(sortedProb, '-r');
hold on;
plot(binCenter, sortedProbBinned, 'xr');
title('New classifier, old features: Probability vs. Frequency in test set');
subplot(2,2,3);
[sortedProb, idx] = sort(test.probRaw);
sortedLabels = test.labels(idx);
sortedProbBinned = arrayfun(@(x,y)sum(sortedLabels(x:y) == 1)./numel(sortedLabels(x:y)), ...
    binLowerLimit, binUpperLimit);
plot(sortedProb, '-g');
hold on;
plot(binCenter, sortedProbBinned, 'xg');
title('New classifier, raw SynEM features: Probability vs. Frequency in test set');
subplot(2,2,4);
[sortedProb, idx] = sort(test.probClass);
sortedLabels = test.labels(idx);
sortedProbBinned = arrayfun(@(x,y)sum(sortedLabels(x:y) == 1)./numel(sortedLabels(x:y)), ...
    binLowerLimit, binUpperLimit);
plot(sortedProb, '-b');
hold on;
plot(binCenter, sortedProbBinned, 'xb');
title('New classifier, raw + class SynEM features: Probability vs. Frequency in test set');

%% Plot precision recall curves
% GP vs. Boosted on old features/rawFeatures/raw+classFeatures

figure;
hold on;

% Old classifier (GP), old features
[Xpr,Ypr,~,AUCpr] = perfcurve(train.labels, train.prob, 1, 'xCrit', 'reca', 'yCrit', 'prec', 'TVals', 0:0.01:1, 'UseNearest', 'off');
plot(Xpr,Ypr, ':k');
label{1} = ['Old classifier (GP), old features, train (AUC: ' num2str(AUCpr) ')'];
[Xpr,Ypr,Tpr,AUCpr] = perfcurve(test.labels, test.prob, 1, 'xCrit', 'reca', 'yCrit', 'prec', 'TVals', 0:0.01:1, 'UseNearest', 'off');
plot(Xpr,Ypr, '-k');
label{2} = ['Old classifier (GP), old features, test (AUC: ' num2str(AUCpr) ')'];
% Additionaly plt 97% precision and accuracy (used for agglomeration in
% last meeting)
plot(Xpr(4), Ypr(4), 'xk');
label{3} = ['Old classifier (GP), old features, test, 97% value used for agglo, prec: ' num2str(Ypr(4)) ', reca: ' num2str(Xpr(4))];

% New classifier (Logit), old features
[Xpr,Ypr,~,AUCpr] = perfcurve(train.labels, train.probOld, 1, 'xCrit', 'reca', 'yCrit', 'prec', 'TVals', 0:0.01:1, 'UseNearest', 'off');
plot(Xpr,Ypr, ':r');
label{4} = ['New classifier (Logit), old features, train (AUC: ' num2str(AUCpr) ')'];
[Xpr,Ypr,~,AUCpr] = perfcurve(test.labels, test.probOld, 1, 'xCrit', 'reca', 'yCrit', 'prec', 'TVals', 0:0.01:1, 'UseNearest', 'off');
plot(Xpr,Ypr, '-r');
label{5} = ['New classifier (Logit), old features, test (AUC: ' num2str(AUCpr) ')'];

% New classifier (Logit), raw SynEM features
[Xpr,Ypr,~,AUCpr] = perfcurve(train.labels, train.probRaw, 1, 'xCrit', 'reca', 'yCrit', 'prec', 'TVals', 0:0.01:1, 'UseNearest', 'off');
plot(Xpr,Ypr, ':g');
label{6} = ['New classifier (Logit), raw SynEM features, train (AUC: ' num2str(AUCpr) ')'];
[Xpr,Ypr,~,AUCpr] = perfcurve(test.labels, test.probRaw, 1, 'xCrit', 'reca', 'yCrit', 'prec', 'TVals', 0:0.01:1, 'UseNearest', 'off');
plot(Xpr,Ypr, '-g');
label{7} = ['New classifier (Logit), raw SynEM features, test (AUC: ' num2str(AUCpr) ')'];

% New classifier (Logit), raw + class SynEM features
[Xpr,Ypr,~,AUCpr] = perfcurve(train.labels, train.probClass, 1, 'xCrit', 'reca', 'yCrit', 'prec', 'TVals', 0:0.01:1, 'UseNearest', 'off');
plot(Xpr,Ypr, ':b');
label{8} = ['New classifier (Logit), raw + class SynEM features, train (AUC: ' num2str(AUCpr) ')'];
[Xpr,Ypr,Tpr,AUCpr] = perfcurve(test.labels, test.probClass, 1, 'xCrit', 'reca', 'yCrit', 'prec', 'TVals', 0:0.01:1, 'UseNearest', 'off');
plot(Xpr,Ypr, '-b');
label{9} = ['New classifier (Logit), raw + class SynEM features, test (AUC: ' num2str(AUCpr) ')'];
% Plot range of 
idx = find(Ypr > 0.99, 1, 'last');
plot(Xpr(idx), Ypr(idx), 'xb');
label{10} = ['New classifier (Logit), raw + class SynEM features, test, upper limit probability to test: ' num2str(Tpr(idx)) ' , prec: ' num2str(Ypr(idx)) ', reca: ' num2str(Xpr(idx))];
idx = find(Ypr > 0.97, 1, 'last');
plot(Xpr(idx), Ypr(idx), 'xb');
label{11} = ['New classifier (Logit), raw + class SynEM features, test, lower limit probability to test: ' num2str(Tpr(idx)) ' , prec: ' num2str(Ypr(idx)) ', reca: ' num2str(Xpr(idx))];


legend(label, 'Location', 'southwest');
xlabel('Recall');
ylabel('Precision');
xlim([0 1]);
ylim([0 1]);
axis square;
