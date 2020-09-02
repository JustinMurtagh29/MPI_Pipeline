% compare feature importance between Mk and ex145 segmentClassifier

rootDir = '/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet';
addpath(genpath('/gaba/u/sahilloo/repos/amotta/matlab/'))
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;
featureSetName = 'segmentAgglomerate';

outDir = fullfile(param.saveFolder,'typeEM','spine', featureSetName);

% ex145 segment
m = load('/gaba/u/mberning/results/pipeline/20170217_ROI/segmentClassifier.mat');
imp = predictorImportance(m.classifiers{4}); % spinehead
sprintf('%d feats out of %d were > 0', sum(imp>0), numel(imp))

figure;
plot(imp,'-or');
hold on;

% MK: trained with typeEM - labels, segAgglo based features
m = load('/u/sahilloo/Mk1_F6_JS_SubI_v1/type-em/spineClassifier/debug.mat');
imp = predictorImportance(m.curClassifier.classifiers{1});
sprintf('%d feats out of %d were > 0', sum(imp>0), numel(imp))

plot(imp, '-+b')
xlabel('Feature index')
ylabel('Feature importance')
legend({'ex145', 'Mk'})
saveas(gcf, fullfile(outDir, 'featImpCompareWithL4.png'))
close all


