%% Settings
% Load parameter of pipeline run
load /gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat;

% Extract information of dense cube annotations
newGTfolder = '+connectEM/trainingData/afterFNandMergerAnnotationMB3/';

% Define bounding box of training regions as passed to annotators
% +1 for wk vs. matlab coordinate system offset
pT.local(1).bboxSmall = [4133 3963 2253; 4578 4408 2432]'+1;
pT.local(2).bboxSmall = [4438 1320 893; 4883 1765 1072]'+1;
pT.local(3).bboxSmall = [1824 6673 1239; 2269 7118 1418]'+1;

pT.local(1).trainFile{1} = [newGTfolder 'region-1.nml'];
pT.local(2).trainFile{1} = [newGTfolder 'region-2.nml'];
pT.local(3).trainFile{1} = [newGTfolder 'region-3.nml'];

segMeta = load([p.saveFolder 'segmentMeta.mat']);
segMeta.point = segMeta.point';

%% Collect training data from segmentation + nmls
% Get a list of labels from each (redundant) training file in each region
gt = connectEM.getContinuityLabelsFromNml(p, pT);

%% Collect features from SynEM feature calculation for all edges in GT
% Note this will only keep borders > 10 voxel and sorts out edges that have
% multiple borders
gt = connectEM.getFeaturesForEdges(p, gt);
gt = rmfield(gt, {'segIdsOfGTnodes', 'segIdsGT', 'mergedSegments', 'leftSegments'});

%% Write FFN and canidates out for annotation (based on edge annotation from now on) again

outputFolder = ['+connectEM' filesep 'trainingData' filesep 'edgeAnnotation1' filesep];
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Now including border node as suggested by Benedikt
for i=1:length(gt)
    idx = find(gt(i).labels == -1 & gt(i).prob > .9);
    nodes = mat2cell([segMeta.point(gt(i).edges(idx,1),:) round(gt(i).borderCoM(idx,:)) segMeta.point(gt(i).edges(idx,2),:)], ...
        ones(size(gt(i).edges(idx,:),1),1), 9);
    nodes = cellfun(@(x)reshape(x,3,3)', nodes, 'uni', 0);
    treeNames = arrayfun(@(x)['region' num2str(i) '_score' num2str(gt(i).prob(x), '%3.2f') '_component' num2str(x, '%.5i') ], idx, 'uni', 0);
    connectEM.generateSkeletonFromNodes([outputFolder 'region' num2str(i) '_fnCanidates.nml'], nodes, treeNames, []);
    idx = find(gt(i).labels == 1 & gt(i).prob < .05);
    nodes = mat2cell([segMeta.point(gt(i).edges(idx,1),:) round(gt(i).borderCoM(idx,:)) segMeta.point(gt(i).edges(idx,2),:)], ...
        ones(size(gt(i).edges(idx,:),1),1), 9);
    nodes = cellfun(@(x)reshape(x,3,3)', nodes, 'uni', 0);
    treeNames = arrayfun(@(x)['region' num2str(i) '_score' num2str(gt(i).prob(x), '%3.2f') '_component' num2str(x, '%.5i') ], idx, 'uni', 0);
    connectEM.generateSkeletonFromNodes([outputFolder 'region' num2str(i) '_fpCanidates.nml'], nodes, treeNames, []);
end
% In case I mess sth. up
%save('/home/mberning/classifierTemp.mat', '-v7.3');

%% Read edge annotation results and incorporate into back into GT
inputFolder = ['+connectEM' filesep 'trainingData' filesep 'edgeAnnotation1result' filesep];

for i=1:length(gt)
    % False positive
    fpSkel = skeleton([inputFolder 'region' num2str(i) '_fpCanidates.nml']);
    fpResult = cellfun(@(x)regexp(x, 'region(\d{1})_score(\d{1}.\d{2})_component(\d{5}).*(\w{2,3})', 'tokens'), fpSkel.names);
    fpResult = cat(1, fpResult{:});
    decision = fpResult(:,4);
    assert(all(strcmp(decision, 'es') | strcmp(decision, 'No')));
    temp = strcmp(decision, 'es');
    decision = zeros(size(temp));
    decision(temp) = 1;
    decision(~temp) = -1;
    assert(all(decision == -1 | decision == 1));
    fpResult = cellfun(@str2double, fpResult(:,1:3));
    assert(all(fpResult(:,1) == i));
    assert(all(fpResult(:,2) == arrayfun(@(x)str2double(num2str(x, '%3.2f')), gt(i).prob(fpResult(:,3)))));
    assert(all(gt(i).labels(fpResult(:,3)) == 1));
    gt(i).labels(fpResult(:,3)) = decision;
    % False negative
    fnSkel = skeleton([inputFolder 'region' num2str(i) '_fnCanidates.nml']);
    fnResult = cellfun(@(x)regexp(x, 'region(\d{1})_score(\d{1}.\d{2})_component(\d{5}).*(\w{2,3})', 'tokens'), fnSkel.names);
    fnResult = cat(1, fnResult{:});
    decision = fnResult(:,4);
    assert(all(strcmp(decision, 'es') | strcmp(decision, 'No')));
    temp = strcmp(decision, 'es');
    decision = zeros(size(temp));
    decision(temp) = 1;
    decision(~temp) = -1;
    assert(all(decision == -1 | decision == 1));
    fnResult = cellfun(@str2double, fnResult(:,1:3));
    assert(all(fnResult(:,1) == i));
    assert(all(fnResult(:,2) == arrayfun(@(x)str2double(num2str(x, '%3.2f')), gt(i).prob(fnResult(:,3)))));
    assert(all(gt(i).labels(fnResult(:,3)) == -1));
    gt(i).labels(fnResult(:,3)) = decision;
    % Display numbers FP and FN remaining:
    display(['Training region' num2str(i)]);
    idx = find(gt(i).labels == 1 & gt(i).prob < .05);
    display(['False positive: ' num2str(numel(idx))]);    
    idx = find(gt(i).labels == -1 & gt(i).prob > .9);
    display(['False negative: ' num2str(numel(idx))]);    
end

%% Write all edges?

%{

%% Write all labels extracted from dense cubes to nml's
for i=1:length(gt)
    posEdges = gt(i).edges(gt(i).labels == 1,:);
    negEdges = gt(i).edges(gt(i).labels == -1,:);
    ccPos = Graph.findConnectedComponents(posEdges, true, true);
    ccNeg = Graph.findConnectedComponents(negEdges, true, true);
    connectEM.generateSkeletonFromAgglo(posEdges, segMeta.point, ccPos, strseq('positiveEdges ', 1:length(ccPos)), ... 
    '+connectEM/trainingData/denseSkel/', segMeta.maxSegId);
    connectEM.generateSkeletonFromAgglo(negEdges, segMeta.point, ccNeg, strseq('negativeEdges ', 1:length(ccNeg)), ...
    '+connectEM/trainingData/denseSkel/', segMeta.maxSegId);
end

%}

%% Divide GT into training an test set (80%/20%) in each dense cube

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

%% Augment training & test set (by adding features for "inverted" edges)

load([p.saveFolder 'SynapseClassifier.mat'], 'fm');
fm.areaT = 10; 
train.classFeaturesInv =  fm.invertDirection(train.classFeatures);
train.rawFeaturesInv =  fm.invertDirection(train.classFeatures);
test.classFeaturesInv =  fm.invertDirection(test.classFeatures);
test.rawFeaturesInv =  fm.invertDirection(test.classFeatures);

%% plot class distribution of some collected features
for i=1:10
    figure('Position', [1 1 1920 999]);
    histogram(train.classFeatures(train.labels == 1,i))
    hold on; histogram(train.classFeatures(train.labels == -1,i))
end

%% Train different classifier on all 3 feature training sets
classifierNewFeatures = connectEM.trainClassifier( ...
    cat(1,cat(2, train.rawFeatures, train.classFeatures),cat(2, train.rawFeaturesInv, train.classFeaturesInv)), ...
    cat(1, train.labels, train.labels));

%% Calculate predicition values for all classifiers (training & test)
[~, train.scoresClass] = classifierNewFeatures.predict( ...
    cat(1,cat(2, train.rawFeatures, train.classFeatures),cat(2, train.rawFeaturesInv, train.classFeaturesInv)));

[~, test.scoresClass] = classifierNewFeatures.predict( ...
    cat(1,cat(2, test.rawFeatures, test.classFeatures),cat(2, test.rawFeaturesInv, test.classFeaturesInv)));

%% Transfer to proabilities
sigmoid = @(x)1./(1+exp(-1.*x));

train.probClass = sigmoid(train.scoresClass(:,1));

test.probClass = sigmoid(test.scoresClass(:,1));

%% augmentation: probability edge vs. inverted/augmented edge
% Load state of training data and predicted proabilities of classifier
% trained afterEdgeAnnotation1 update for one more round of edge updates
load('/gaba/u/mberning/results/edgeClassifier/stateAfterEdgeAnnotation1.mat');

lastElement = size(test.probClass,1)./2;
plot(test.probClass(1:lastElement), test.probClass(lastElement+1:end), 'xb');
title('Test set: Edge vs. inverted edge');
xlabel('Probability edge');
ylabel('Probability inverted edge');

lastElement = size(train.probClass,1)./2;
plot(train.probClass(1:lastElement), train.probClass(lastElement+1:end), 'xb');
title('Training set: Edge vs. inverted edge');
xlabel('Probability edge');
ylabel('Probability inverted edge');

%% comparison to last classifier

lastElement = size(test.probClass,1)./2;
plot(test.probClass(1:lastElement), test.prob, 'xb');
title('Test set: Changes of probabilities wrt to old classifier');
xlabel('Probability edge');
ylabel('Probability edge last classifier');

lastElement = size(train.probClass,1)./2;
plot(train.probClass(1:lastElement), train.prob, 'xb');
title('Training set: Changes of probabilities wrt to old classifier');
xlabel('Probability edge');
ylabel('Probability edge last classifier');

%% Write FN and FP for annotation (based on edge annotation from now on) again

outputFolder = ['+connectEM' filesep 'trainingData' filesep 'edgeAnnotation2' filesep];
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Now including border node as suggested by Benedikt (and connecting nodes
% linearly now)

% test
i = 1;
idx = find(test(i).labels == -1 & test(i).probClass(1:length(test(i).probClass)/2) > .9);
nodes = mat2cell([segMeta.point(test(i).edges(idx,1),:) round(test(i).borderCoM(idx,:)) segMeta.point(test(i).edges(idx,2),:)], ...
    ones(size(test(i).edges(idx,:),1),1), 9);
nodes = cellfun(@(x)reshape(x,3,3)', nodes, 'uni', 0);
treeNames = arrayfun(@(x)['test_score' num2str(test(i).probClass(x), '%3.2f') '_component' num2str(x, '%.5i') ], idx, 'uni', 0);
connectEM.generateSkeletonFromNodes([outputFolder 'test_fnCanidates.nml'], nodes, treeNames, [], true);
idx = find(test(i).labels == 1 & test(i).probClass(1:length(test(i).probClass)/2) < .1);
nodes = mat2cell([segMeta.point(test(i).edges(idx,1),:) round(test(i).borderCoM(idx,:)) segMeta.point(test(i).edges(idx,2),:)], ...
    ones(size(test(i).edges(idx,:),1),1), 9);
nodes = cellfun(@(x)reshape(x,3,3)', nodes, 'uni', 0);
treeNames = arrayfun(@(x)['test_score' num2str(test(i).probClass(x), '%3.2f') '_component' num2str(x, '%.5i') ], idx, 'uni', 0);
connectEM.generateSkeletonFromNodes([outputFolder 'test_fpCanidates.nml'], nodes, treeNames, [], true);
% training
i = 1;
idx = find(train(i).labels == -1 & train(i).probClass(1:length(train(i).probClass)/2) > .9);
nodes = mat2cell([segMeta.point(train(i).edges(idx,1),:) round(train(i).borderCoM(idx,:)) segMeta.point(train(i).edges(idx,2),:)], ...
    ones(size(train(i).edges(idx,:),1),1), 9);
nodes = cellfun(@(x)reshape(x,3,3)', nodes, 'uni', 0);
treeNames = arrayfun(@(x)['train_score' num2str(train(i).probClass(x), '%3.2f') '_component' num2str(x, '%.5i') ], idx, 'uni', 0);
connectEM.generateSkeletonFromNodes([outputFolder 'train_fnCanidates.nml'], nodes, treeNames, [], true);
idx = find(train(i).labels == 1 & train(i).probClass(1:length(train(i).probClass)/2) < .1);
nodes = mat2cell([segMeta.point(train(i).edges(idx,1),:) round(train(i).borderCoM(idx,:)) segMeta.point(train(i).edges(idx,2),:)], ...
    ones(size(train(i).edges(idx,:),1),1), 9);
nodes = cellfun(@(x)reshape(x,3,3)', nodes, 'uni', 0);
treeNames = arrayfun(@(x)['train_score' num2str(train(i).probClass(x), '%3.2f') '_component' num2str(x, '%.5i') ], idx, 'uni', 0);
connectEM.generateSkeletonFromNodes([outputFolder 'train_fpCanidates.nml'], nodes, treeNames, [], true);

%% most changed probabilities
idx = find(abs(test(i).probClass(1:length(test(i).probClass)/2) - test(i).prob) > 0.5);
nodes = mat2cell([segMeta.point(test(i).edges(idx,1),:) round(test(i).borderCoM(idx,:)) segMeta.point(test(i).edges(idx,2),:)], ...
    ones(size(test(i).edges(idx,:),1),1), 9);
nodes = cellfun(@(x)reshape(x,3,3)', nodes, 'uni', 0);
treeNames = arrayfun(@(x)['region' num2str(i) '_score' num2str(test(i).probClass(x), '%3.2f') '_scoreOld' num2str(test(i).prob(x), '%3.2f') ], idx, 'uni', 0);
connectEM.generateSkeletonFromNodes([outputFolder 'region' num2str(i) 'changedProbs.nml'], nodes, treeNames, []);


%% Visualize predicted probabilities vs. frequencies in test set
binSize = 500;
binLowerLimit = 1:binSize:length(test.probClass)/2;
binUpperLimit = [(binSize):binSize:length(test.probClass)/2 length(test.probClass)/2];
binCenter = (binLowerLimit + binUpperLimit) ./ 2;
figure;
% subplot(2,2,1);
% [sortedProb, idx] = sort(test.prob);
% sortedLabels = test.labels(idx);
% sortedProbBinned = arrayfun(@(x,y)sum(sortedLabels(x:y) == 1)./numel(sortedLabels(x:y)), ...
%     binLowerLimit, binUpperLimit);
% plot(sortedProb, '-k');
% hold on;
% plot(binCenter, sortedProbBinned, 'xk');
% title('Old classifier, old features: Probability vs. Frequency in test set');
% subplot(2,2,2);
% [sortedProb, idx] = sort(test.probOld);
% sortedLabels = test.labels(idx);
% sortedProbBinned = arrayfun(@(x,y)sum(sortedLabels(x:y) == 1)./numel(sortedLabels(x:y)), ...
%     binLowerLimit, binUpperLimit);
% plot(sortedProb, '-r');
% hold on;
% plot(binCenter, sortedProbBinned, 'xr');
% title('New classifier, old features: Probability vs. Frequency in test set');
% subplot(2,2,3);
% [sortedProb, idx] = sort(test.probRaw);
% sortedLabels = test.labels(idx);
% sortedProbBinned = arrayfun(@(x,y)sum(sortedLabels(x:y) == 1)./numel(sortedLabels(x:y)), ...
%     binLowerLimit, binUpperLimit);
% plot(sortedProb, '-g');
% hold on;
% plot(binCenter, sortedProbBinned, 'xg');
% title('New classifier, raw SynEM features: Probability vs. Frequency in test set');
% subplot(2,2,4);
[sortedProb, idx] = sort(test.probClass(1:length(test.probClass)/2));
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

% % Old classifier (GP), old features
% [Xpr,Ypr,~,AUCpr] = perfcurve(train.labels, train.prob, 1, 'xCrit', 'reca', 'yCrit', 'prec', 'TVals', 0:0.001:1);
% Ypr(1) = 1;
% plot(Xpr,Ypr, ':k');
% label{1} = ['Old classifier (GP), old features, train (AUC: ' num2str(AUCpr) ')'];
% [Xpr,Ypr,Tpr,AUCpr] = perfcurve(test.labels, test.prob, 1, 'xCrit', 'reca', 'yCrit', 'prec', 'TVals', 0:0.001:1);
% Ypr(1) = 1;
% plot(Xpr,Ypr, '-k');
% label{2} = ['Old classifier (GP), old features, test (AUC: ' num2str(AUCpr) ')'];
% % Additionaly plt 97% precision and accuracy (used for agglomeration in
% % last meeting)
% [~, idx] = min(abs(Tpr - 0.97));
% plot(Xpr(idx), Ypr(idx), 'xk');
% label{3} = ['Old classifier (GP), old features, test, 97% value used for agglo, prec: ' num2str(Ypr(idx)) ', reca: ' num2str(Xpr(idx))];

% % New classifier (Logit), old features
% [Xpr,Ypr,~,AUCpr] = perfcurve(train.labels, train.probOld, 1, 'xCrit', 'reca', 'yCrit', 'prec', 'TVals', 0:0.001:1);
% Ypr(1) = 1;
% plot(Xpr,Ypr, ':r');
% label{4} = ['New classifier (Logit), old features, train (AUC: ' num2str(AUCpr) ')'];
% [Xpr,Ypr,~,AUCpr] = perfcurve(test.labels, test.probOld, 1, 'xCrit', 'reca', 'yCrit', 'prec', 'TVals', 0:0.001:1);
% Ypr(1) = 1;
% plot(Xpr,Ypr, '-r');
% label{5} = ['New classifier (Logit), old features, test (AUC: ' num2str(AUCpr) ')'];
% 
% % New classifier (Logit), raw SynEM features
% [Xpr,Ypr,~,AUCpr] = perfcurve(train.labels, train.probRaw, 1, 'xCrit', 'reca', 'yCrit', 'prec', 'TVals', 0:0.001:1);
% Ypr(1) = 1;
% plot(Xpr,Ypr, ':g');
% label{6} = ['New classifier (Logit), raw SynEM features, train (AUC: ' num2str(AUCpr) ')'];
% [Xpr,Ypr,~,AUCpr] = perfcurve(test.labels, test.probRaw, 1, 'xCrit', 'reca', 'yCrit', 'prec', 'TVals', 0:0.001:1);
% Ypr(1) = 1;
% plot(Xpr,Ypr, '-g');
% label{7} = ['New classifier (Logit), raw SynEM features, test (AUC: ' num2str(AUCpr) ')'];

% New classifier (Logit), raw + class SynEM features
[Xpr,Ypr,~,AUCpr] = perfcurve(cat(1,train.labels,train.labels), train.probClass, 1, 'xCrit', 'reca', 'yCrit', 'prec', 'TVals', 0:0.001:1);
Ypr(1) = 1;
plot(Xpr,Ypr, ':b');
label{1} = ['New classifier (Logit), raw + class SynEM features, train (AUC: ' num2str(AUCpr) ')'];
[Xpr,Ypr,Tpr,AUCpr] = perfcurve(cat(1, test.labels, test.labels), test.probClass, 1, 'xCrit', 'reca', 'yCrit', 'prec', 'TVals', 0:0.001:1);
Ypr(1) = 1;
plot(Xpr,Ypr, '-b');
label{2} = ['New classifier (Logit), raw + class SynEM features, test (AUC: ' num2str(AUCpr) ')'];
% Plot range of 
idx = find(Ypr > 0.95, 1, 'last');
plot(Xpr(idx), Ypr(idx), 'xb');
label{3} = ['New classifier (Logit), raw + class SynEM features, test, upper limit probability to test: ' num2str(Tpr(idx)) ' , prec: ' num2str(Ypr(idx)) ', reca: ' num2str(Xpr(idx))];
idx = find(Ypr > 0.94, 1, 'last');
plot(Xpr(idx), Ypr(idx), 'xb');
label{4} = ['New classifier (Logit), raw + class SynEM features, test, lower limit probability to test: ' num2str(Tpr(idx)) ' , prec: ' num2str(Ypr(idx)) ', reca: ' num2str(Xpr(idx))];


legend(label, 'Location', 'southwest');
xlabel('Recall');
ylabel('Precision');
xlim([0 1]);
ylim([0 1]);
axis square;

%% Determine whether compacted classifier works here
classifier = compact(classifierNewFeatures);
save(['/gaba/u/mberning/results/edgeClassifier/' datestr(clock,30) '.mat'], 'classifier');

a = classifier.predict(cat(2, test.rawFeatures, ...
test.classFeatures));
b = classifierNewFeatures.predict(cat(2, test.rawFeatures, ...
test.classFeatures));

all(a(:) == b(:))

%%  Comparison of performance old classifier on new test set
a = load('/gaba/u/mberning/results/edgeClassifier/20170210T121156.mat', 'classifier');

[~, test.scoresClassOld] = a.classifier.predict( ...
    cat(2, test.rawFeatures, test.classFeatures));
[~, train.scoresClassOld] = a.classifier.predict( ...
    cat(2, train.rawFeatures, train.classFeatures));

test.probClassOld = sigmoid(test.scoresClassOld(:,1));
train.probClassOld = sigmoid(train.scoresClassOld(:,1));

%% Plotting
figure;
hold on;

[Xpr,Ypr,~,AUCpr] = perfcurve(test.labels, test.probClassOld, 1, 'xCrit', 'reca', 'yCrit', 'prec', 'TVals', 0:0.001:1);
Ypr(1) = 1;
plot(Xpr,Ypr, '-g');
label{1} = ['Old classifier (the one with the good PR curve before) on new GT, test (AUC: ' num2str(AUCpr) ')'];
[Xpr,Ypr,Tpr,AUCpr] = perfcurve(cat(1, test.labels, test.labels), test.probClass, 1, 'xCrit', 'reca', 'yCrit', 'prec', 'TVals', 0:0.001:1);
Ypr(1) = 1;
plot(Xpr,Ypr, '-b');
label{2} = ['New classifier on new GT, test (AUC: ' num2str(AUCpr) ')'];
% Plot range of 
idx = find(Ypr > 0.95, 1, 'last');
plot(Xpr(idx), Ypr(idx), 'xb');
label{3} = ['New classifier (Logit), raw + class SynEM features, test, upper limit probability to test: ' num2str(Tpr(idx)) ' , prec: ' num2str(Ypr(idx)) ', reca: ' num2str(Xpr(idx))];
idx = find(Ypr > 0.94, 1, 'last');
plot(Xpr(idx), Ypr(idx), 'xb');
label{4} = ['New classifier (Logit), raw + class SynEM features, test, lower limit probability to test: ' num2str(Tpr(idx)) ' , prec: ' num2str(Ypr(idx)) ', reca: ' num2str(Xpr(idx))];

legend(label, 'Location', 'southwest');
xlabel('Recall');
ylabel('Precision');
xlim([0 1]);
ylim([0 1]);
axis square;

%%  Comparison of performance old classifier on old test set
load /home/mberning/Desktop/oldGT.mat;

figure;
hold on;

labels = cat(2, gtOld(:).labels)';
[oldEdges, idx] = unique(cat(1, gtOld(:).edges), 'rows');
labels = labels(idx);

prob = cat(1, train.probClass(1:size(train.probClass)/2), test.probClass(1:size(test.probClass)/2));
probOld = cat(1, train.probClassOld, test.probClassOld);
[newEdges, idx] = unique(cat(1, train.edges, test.edges), 'rows');
prob = prob(idx);
probOld = probOld(idx);
edges = intersect(newEdges, oldEdges, 'rows');
idx1 = ismember(oldEdges, edges, 'rows');
idx2 = ismember(newEdges, edges, 'rows');
labels = labels(idx1);
prob = prob(idx2);
probOld = probOld(idx2);

[Xpr,Ypr,~,AUCpr] = perfcurve(labels, probOld, 1, 'xCrit', 'reca', 'yCrit', 'prec', 'TVals', 0:0.001:1);
Ypr(1) = 1;
plot(Xpr,Ypr, '-g');
label{1} = ['Old classifier (the one with the good PR curve before) on old GT, test (AUC: ' num2str(AUCpr) ')'];
[Xpr,Ypr,Tpr,AUCpr] = perfcurve(labels, prob, 1, 'xCrit', 'reca', 'yCrit', 'prec', 'TVals', 0:0.001:1);
Ypr(1) = 1;
plot(Xpr,Ypr, '-b');
label{2} = ['New classifier on old GT, test (AUC: ' num2str(AUCpr) ')'];
% Plot range of 
idx = find(Ypr > 0.95, 1, 'last');
plot(Xpr(idx), Ypr(idx), 'xb');
label{3} = ['New classifier (Logit), raw + class SynEM features, test, upper limit probability to test: ' num2str(Tpr(idx)) ' , prec: ' num2str(Ypr(idx)) ', reca: ' num2str(Xpr(idx))];
idx = find(Ypr > 0.94, 1, 'last');
plot(Xpr(idx), Ypr(idx), 'xb');
label{4} = ['New classifier (Logit), raw + class SynEM features, test, lower limit probability to test: ' num2str(Tpr(idx)) ' , prec: ' num2str(Ypr(idx)) ', reca: ' num2str(Xpr(idx))];

legend(label, 'Location', 'southwest');
xlabel('Recall');
ylabel('Precision');
xlim([0 1]);
ylim([0 1]);
axis square;

%% Debug something

plot(train.scoresClassOld, train.scoresClass(1:size(train.scoresClass,1)/2), 'x');
plot(test.scoresClassOld, test.scoresClass(1:size(test.scoresClass,1)/2), 'x');
plot(train.probClassOld, train.probClass(1:size(train.probClass,1)/2), 'x');
plot(test.probClassOld, test.probClass(1:size(test.probClass,1)/2), 'x');
