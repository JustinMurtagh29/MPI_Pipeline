function y = evalutateAggloMetaMeta(graph, agglos)
addpath(genpath('/gaba/u/kboerg/code/auxiliaryMethods'))
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameter.mat', 'p')
segmentMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta.mat')
evalfolder = '/gaba/u/kboerg/code/pipeline/+connectEM/evaluationData/'
skelpath = {[evalfolder '2012-09-28_ex145_07x2_new2__explorational__mhelmstaedter__1bad47-2.nml']};
meta_liste = {'dendrite_gt', 'new_axon_gt', 'old_axon_gt'};
for meta_idx = 1 : 3
    liste = dir([evalfolder meta_liste{meta_idx} '/*.nml']);
    for file_idx = 1 : length(liste)
        skelpath{end + 1} = [evalfolder meta_liste{meta_idx} filesep liste(file_idx).name];
    end
end
y = connectEM.evaluateAggloMeta(skelpath, graph, segmentMeta, agglos, p)
