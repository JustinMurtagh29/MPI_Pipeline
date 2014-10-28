function classifierVisualization(groundTruthFile, hyperFile, outputDir)

% load data
load(groundTruthFile); % training data generated with prepareTrainingData
load(hyperFile); % all parameter for GP (e.g. hyp, meanfunc)

% gpml toolbox usage
run('/zdata/manuel/code/active/gpml/startup.m');
% Make predictions
[trMean trVar trLatMean trLatVar trLp] = gp(hyp, inffunc, meanfunc, covfunc, likfunc, trainingData, trainingLabels, trainingData, ones(size(trainingData,1), 1));
[teMean teVar teLatMean teLatVar teLp] = gp(hyp, inffunc, meanfunc, covfunc, likfunc, trainingData, trainingLabels, testData, ones(size(testData,1), 1));
trProb = exp(trLp);
teProb = exp(teLp);

% Save data in case anything goes wrong, calculation above takes quite some time
save('/zdata/manuel/sync/classifier/data.mat');

% Visualization, each plot inclunding training and test data to judge overrfitting/model problems

% Plot histogramms of predictions for different classes
figure('Units', 'centimeters', 'Position', [0 0 29.7 21], 'Visible', 'off', ...
    'PaperType', 'A4', 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 29.7 21], 'PaperOrientation', 'portrait');
a = subplot(2,2,1);
hist(teProb(testLabels == 1),0.005:0.01:0.995);
title('Probabilities for positive test labels');
b = subplot(2,2,2);
hist(teProb(testLabels == -1),0.005:0.01:0.995);
title('Probabilities for negative test labels');
c = subplot(2,2,3);
hist(trProb(trainingLabels == 1),0.005:0.01:0.995);
title('Probabilities for positive test labels');
d = subplot(2,2,4);
hist(trProb(trainingLabels == -1),0.005:0.01:0.995);
title('Probabilities for negative test labels');
linkaxes([a b c d], 'xy');
set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize')));
print(gcf, [outputDir 'histograms.pdf', '-r600', '-dpdf');

% Plot boxplot
figure('Units', 'centimeters', 'Position', [0 0 29.7 21], 'Visible', 'off', ...
    'PaperType', 'A4', 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 29.7 21], 'PaperOrientation', 'portrait');
a = subplot(2,1,1);
boxplot(teProb,testLabels);
title('Predicted probabilities for both classes in test data');
b = subplot(2,1,2);
boxplot(trProb,trainingLabels);
title('Predicted probabilities for both classes in training data');
linkaxes([a b], 'xy');
set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize')));
print(gcf, [outputDir 'boxplot.pdf'], '-r600', '-dpdf');

% Plot Precission-Recall
for i=1;100

figure('Units', 'centimeters', 'Position', [0 0 29.7 21], 'Visible', 'off', ...
    'PaperType', 'A4', 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 29.7 21], 'PaperOrientation', 'portrait');
a = subplot(2,1,1);
boxplot(teProb,testLabels);
title('Predicted probabilities for both classes in test data');
b = subplot(2,1,2);
boxplot(trProb,trainingLabels);
title('Predicted probabilities for both classes in training data');
linkaxes([a b], 'xy');
set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize')));
print(gcf, [outputDir 'precisionRecall.pdf'], '-r600', '-dpdf');

end

