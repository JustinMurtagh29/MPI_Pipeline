% compare feature importance between Mk and ex145 segmentClassifier

rootDir = '/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet';
addpath(genpath('/gaba/u/sahilloo/repos/amotta/matlab/'))
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;
featureSetName = 'segmentAgglomerate';

outDir = fullfile(param.saveFolder,'typeEM','spine', featureSetName);

% ex145 segment
m = load('/gaba/u/mberning/results/pipeline/20170217_ROI/segmentClassifier.mat');
imp1 = predictorImportance(m.classifiers{4}); % spinehead
[imp1,idxSort] = sort(imp1);

sprintf('%d feats out of %d were > 0', sum(imp1>0), numel(imp1))
sprintf('Mean of non-zero: %f', mean(nonzeros(imp1)))
sprintf('Median of non-zero: %f', median(nonzeros(imp1)))

figure;
bar(imp1,'EdgeColor','r', 'FaceColor', 'none');
hold on;

% MK: trained with typeEM - labels, segAgglo based features
m = load('/u/sahilloo/Mk1_F6_JS_SubI_v1/type-em/spineClassifier/debug.mat');
imp2 = predictorImportance(m.curClassifier.classifiers{1});
imp2 = imp2(idxSort);

sprintf('%d feats out of %d were > 0', sum(imp2>0), numel(imp2))
sprintf('Mean of non-zero: %f', mean(nonzeros(imp2)))
sprintf('Median of non-zero: %f', median(nonzeros(imp2)))

bar(imp2, 'EdgeColor', 'b', 'FaceColor', 'none')
xlabel('Feature index')
ylabel('Feature importance')
legend({'ex145', 'Mk'})
saveas(gcf, fullfile(outDir, 'featImpCompareWithL4.png'))
close all


