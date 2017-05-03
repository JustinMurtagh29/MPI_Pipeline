function y = evaluateAggloMetaMeta(graph, agglos_axon, agglos_dendrite, nameOfAgglo, segmentMeta)

% Load pipeline information
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameter.mat', 'p');
% Where GT can be found (not sure why called "evalfolder")
evalfolder = '/gaba/u/kboerg/code/pipeline/+connectEM/evaluationData/';
% File locations of GT
skelpath = {[evalfolder '2012-09-28_ex145_07x2_new2__explorational__mhelmstaedter__1bad47-2_axon.nml'], ...
    [evalfolder '2012-09-28_ex145_07x2_ROI2017__explorational__amotta__b8e060.nml'], ...
    [evalfolder '2012-09-28_ex145_07x2_new2__explorational__mhelmstaedter__b8e053_dend1.nml'], ...
    [evalfolder '2012-09-28_ex145_07x2_new2__explorational__mhelmstaedter__b8e053_dend2.nml'], ...
    [evalfolder '2012-09-28_ex145_07x2_new2__explorational__mhelmstaedter__1bad47-2_dend.nml']};

% Additional GT in subfolders
meta_liste = {'dendrite_gt', 'new_axon_gt', 'old_axon_gt'};
for meta_idx = 1 : 3
    liste = dir([evalfolder meta_liste{meta_idx} '/*.nml']);
    for file_idx = 1 : length(liste)
        skelpath{end + 1} = [evalfolder meta_liste{meta_idx} filesep liste(file_idx).name];
    end
end
skelpath = [skelpath, {[evalfolder '2012-09-28_ex145_07x2_new2__explorational__kboergens__bd3e11.nml'], ...
    [evalfolder '2012-09-28_ex145_07x2_new2__explorational__kboergens__bde7d2.nml'], ...
    [evalfolder '2012-09-28_ex145_07x2_new2__explorational__kboergens__be0e2a.nml'], ...
    [evalfolder '2012-09-28_ex145_07x2_new2__explorational__kboergens__be5275.nml'], ...
    [evalfolder '2012-09-28_ex145_07x2_ROI2017__explorational__amotta__fe80a9.nml'], ...
    [evalfolder '2012-09-28_ex145_07x2_ROI2017__explorational__amotta__f70ad6_axon1.nml'], ...
    [evalfolder '2012-09-28_ex145_07x2_ROI2017__explorational__amotta__f70ad6_axon2.nml'], ...
    [evalfolder '2012-09-28_ex145_07x2_new2__explorational__mhelmstaedter__b8e053_axon.nml']}];



% Where to store results
mainFolder = ['/gaba/scratch/kboerg/eval_agglo/' nameOfAgglo '/'];
mkdir(mainFolder);
% Choose which GT in skelpath to use for evaluation
dendrite_selector = [2:5, 36:39];
axon_selector = 16:25;%[1, 26, 27, 34, 35, 40:43];
% Dendrite
maxTubeDendrite = 5000;
% agglos_dendrite_reverse = createLookup(segmentMeta, agglos_dendrite);
% y.dendrite1 = connectEM.evaluateAggloMeta(skelpath(dendrite_selector), graph, segmentMeta, agglos_dendrite, p, 'dendrites1', mainFolder, 0, 2.5E4, agglos_dendrite_reverse, maxTubeDendrite);
% y.dendrite4 = connectEM.evaluateAggloMeta(skelpath(dendrite_selector), graph, segmentMeta, agglos_dendrite, p, 'dendrites4', mainFolder, 3, 2.5E4, agglos_dendrite_reverse, maxTubeDendrite);
% y.dendritePercolators = cellfun(@(x)sum(segmentMeta.voxelCount(x)), agglos_dendrite(1:min(10,numel(agglos_dendrite))));
% Axon evaluation
maxTubeAxon = 3000;
agglos_axon_reverse = createLookup(segmentMeta, agglos_axon);
y.axon1 = connectEM.evaluateAggloMeta(skelpath(axon_selector), graph, segmentMeta, agglos_axon, p, 'axons1', mainFolder, 0, 6250, agglos_axon_reverse,maxTubeAxon);
%y.axon2 = connectEM.evaluateAggloMeta(skelpath(axon_selector), graph, segmentMeta, agglos_axon, p, 'axons2', mainFolder, 1, 6250, agglos_axon_reverse,maxTubeAxon);
y.axonPercolators = cellfun(@(x)sum(segmentMeta.voxelCount(x)), agglos_axon(1:min(10,numel(agglos_axon))));

end

function agglos_reverse = createLookup(segmentMeta, agglos)

agglos_reverse = zeros(size(segmentMeta.point, 1), 1);
for idx = 1 : length(agglos)
    agglos_reverse(agglos{idx}) = idx;
end

end
