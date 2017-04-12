function y = evalutateAggloMetaMeta(graph, agglos_axon, agglos_dendrite)
addpath(genpath('/gaba/u/kboerg/code/auxiliaryMethods'))
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameter.mat', 'p')
segmentMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta.mat')
evalfolder = '/gaba/u/kboerg/code/pipeline/+connectEM/evaluationData/'
skelpath = {[evalfolder '2012-09-28_ex145_07x2_new2__explorational__mhelmstaedter__1bad47-2_dend.nml'], ...
[evalfolder '2012-09-28_ex145_07x2_new2__explorational__mhelmstaedter__1bad47-2_axon.nml']};
meta_liste = {'dendrite_gt', 'new_axon_gt', 'old_axon_gt'};
for meta_idx = 1 : 3
    liste = dir([evalfolder meta_liste{meta_idx} '/*.nml']);
    for file_idx = 1 : length(liste)
        skelpath{end + 1} = [evalfolder meta_liste{meta_idx} filesep liste(file_idx).name];
    end
end
dendrite_selector = [1, 3:12];
axon_selector = setdiff(1:32, dendrite_selector);

y.axon = connectEM.evaluateAggloMeta(skelpath(axon_selector), graph, segmentMeta, agglos_axon, p)
y.dendrite = connectEM.evaluateAggloMeta(skelpath(dendrite_selector), graph, segmentMeta, agglos_dendrite, p)
