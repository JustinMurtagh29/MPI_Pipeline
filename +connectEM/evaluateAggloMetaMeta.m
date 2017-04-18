function y = evalutateAggloMetaMeta(graph, agglos_axon, agglos_dendrite, nameOfAgglo, segmentMeta)
addpath(genpath('/gaba/u/kboerg/code/auxiliaryMethods'))
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameter.mat', 'p')
%segmentMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta.mat')
evalfolder = '/gaba/u/kboerg/code/pipeline/+connectEM/evaluationData/'
skelpath = {[evalfolder '2012-09-28_ex145_07x2_new2__explorational__mhelmstaedter__1bad47-2_axon.nml'], ...
[evalfolder '2012-09-28_ex145_07x2_ROI2017__explorational__amotta__b8e060.nml'], ...
[evalfolder '2012-09-28_ex145_07x2_new2__explorational__mhelmstaedter__b8e053_dend1.nml'], ...
[evalfolder '2012-09-28_ex145_07x2_new2__explorational__mhelmstaedter__b8e053_dend2.nml'], ...
[evalfolder '2012-09-28_ex145_07x2_new2__explorational__mhelmstaedter__1bad47-2_dend.nml']};
meta_liste = {'dendrite_gt', 'new_axon_gt', 'old_axon_gt'};
for meta_idx = 1 : 3
    liste = dir([evalfolder meta_liste{meta_idx} '/*.nml']);
    for file_idx = 1 : length(liste)
        skelpath{end + 1} = [evalfolder meta_liste{meta_idx} filesep liste(file_idx).name];
    end
end
mainFolder = ['/gaba/scratch/kboerg/eval_agglo/' nameOfAgglo '/'];
mkdir(mainFolder);
dendrite_selector = [2 : 5];
axon_selector = [1 16:35];
%y.dendrite1 = connectEM.evaluateAggloMeta(skelpath(dendrite_selector), graph, segmentMeta, agglos_dendrite, p, 'dendrites1', mainFolder, 0, 1E4);
%y.dendrite4 = connectEM.evaluateAggloMeta(skelpath(dendrite_selector), graph, segmentMeta, agglos_dendrite, p, 'dendrites4', mainFolder, 3, 1E4);
%y.dendrite.percolators = cellfun(@(x)sum(segmentMeta.voxelCount(x)), agglos_dendrite(1:10));
display('reverting agglo')
agglos_axon_reverse = zeros(1, size(segmentMeta.point, 2));
for idx = 1 : length(agglos_axon)
    agglos_axon_reverse(agglos_axon{idx}) = idx;
end

y.axon = connectEM.evaluateAggloMeta(skelpath(axon_selector), graph, segmentMeta, agglos_axon, p, 'axons', mainFolder, 0, 0, agglos_axon_reverse)
