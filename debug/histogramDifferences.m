pAnjali = load('/gaba/u/ganja/results_P14_L4_Corrected/allParameter.mat');
pKevin = load('/gaba/u/mberning/results/pipeline/20151111T183414/allParameter.mat');

normA = load('/gaba/u/sahilloo/gitlab/results/state/normValuesAnjali.mat');
normK = load(pKevin.p.gp.normValues);

weightsA = load(pAnjali.p.local(4,1,4).weightFile);
weightsK = load(pKevin.p.local(4,1,4).weightFile);

weightsAN = normalizeDataForGP(weightsA.weights, false, '/gaba/u/sahilloo/gitlab/results/state/normValuesAnjali.mat');
weightsKN = normalizeDataForGP(weightsK.weights, false, pKevin.p.gp.normValues);

nrFeature = 1;
x = linspace(min([weightsA.weights(:,nrFeature); weightsK.weights(:,nrFeature); weightsAN(:,nrFeature); weightsKN(:,nrFeature)]),...
    max([weightsA.weights(:,nrFeature); weightsK.weights(:,nrFeature); weightsAN(:,nrFeature); weightsKN(:,nrFeature)]),100);
a = hist(weightsA.weights(:,nrFeature), x) ./ numel(weightsA.weights(:,nrFeature));
b = hist(weightsK.weights(:,nrFeature), x) ./ numel(weightsK.weights(:,nrFeature)); 
c = hist(weightsAN(:,nrFeature), x) ./ numel(weightsAN(:,nrFeature));
d = hist(weightsKN(:,nrFeature), x) ./ numel(weightsKN(:,nrFeature)); 
bar([a; b; c; d]'); title(weightsA.weightNames{nrFeature}); legend({'Anjali' 'Kevin' 'Anjali norm' 'Kevin norm'});
set(gca, 'XTick', 1:10:length(x));
set(gca, 'XTickLabel', arrayfun(@num2str, x(1:10:end), 'uni', 0));

