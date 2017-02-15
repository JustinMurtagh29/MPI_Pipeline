diary(['/home/mberning/Desktop/newContinuity/log' datestr(clock, 30) '.txt']);

display(['----- Old training data -----']);

load /gaba/u/mberning/results/pipeline/20161120_ROI/allParameter.mat;

pT.local(1).bboxSmall = [4133 3963 2253; 4578 4408 2432]';
pT.local(2).bboxSmall = [4438 1320 893; 4883 1765 1072]';
pT.local(3).bboxSmall = [1824 6673 1239; 2269 7118 1418]';

pT.local(1).trainFile = {'/home/mberning/Downloads/region1NEW.nml'};
pT.local(2).trainFile = {'/home/mberning/Downloads/region2NEW.nml'};
pT.local(3).trainFile = {'/home/mberning/Downloads/region3NEW.nml'};

gt = connectEM.getContinuityLabelsFromNml(p, pT);

%% Visualize precision recall of old interface classifier

labels = cat(1, gt(:).labels);
prob = cat(1, gt(:).prob);

figure;
[Xpr,Ypr,Tpr,AUCpr] = perfcurve(labels(labels ~= 0), prob(labels ~= 0), 1, 'xCrit', 'reca', 'yCrit', 'prec');
plot(Xpr,Ypr);
xlabel('Recall'); ylabel('Precision');
xlim([0 1]); ylim([0 1]);
title(['Precision-recall curve (AUC: ' num2str(AUCpr) ')']);
saveas(gcf, '/home/mberning/Desktop/newContinuity/prOldGTOldClass.png');

%% New training data

display(['----- New training data -----']);

newGTfolder = '+connectEM/trainingData/';
files = dir([newGTfolder '*.nml']);
temp = cellfun(@(x)strsplit(x, '__'), {files(:).name}, 'uni', 0);
taskIds = cellfun(@(x)x{2}, temp, 'uni', 0);
tracer = cellfun(@(x)x{3}, temp, 'uni', 0);
annotationIds = cellfun(@(x)x{4}(1:end-4), temp, 'uni', 0);

files = {files(:).name};
% Note has to fit bboxSmall, in order only by luck
pT.local(1).trainFile = cellfun(@(x)[newGTfolder x], files(1:3), 'uni', 0);
pT.local(2).trainFile = cellfun(@(x)[newGTfolder x], files(4:6), 'uni', 0);
pT.local(3).trainFile = cellfun(@(x)[newGTfolder x], files(7:9), 'uni', 0);

gtNew = connectEM.getContinuityLabelsFromNml(p, pT, false);

%% Old classifier on new training data

labels = cat(1, gtNew(:).labels);
prob = cat(1, gtNew(:).prob);

figure;
[Xpr,Ypr,Tpr,AUCpr] = perfcurve(labels(labels ~= 0), prob(labels ~= 0), 1, 'xCrit', 'reca', 'yCrit', 'prec');
plot(Xpr,Ypr);
xlabel('Recall'); ylabel('Precision');
xlim([0 1]); ylim([0 1]);
title(['Precision-recall curve (AUC: ' num2str(AUCpr) ')']);
saveas(gcf, '/home/mberning/Desktop/newContinuity/prNewGTOldClass.png');

diary off;