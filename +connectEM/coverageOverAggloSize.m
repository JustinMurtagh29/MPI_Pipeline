
% GT for evalulation
% Where GT can be found (not sure why called "evalfolder")
evalfolder = '/gaba/u/kboerg/code/pipeline/+connectEM/evaluationData/';
mainFolder = '/gaba/scratch/mberning/eval_agglo/';
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
axon_selector = 16:25;

% 
agglos_axon = agglos;
agglos_axon_reverse = zeros(size(segmentMeta.point, 1), 1);
for idx = 1 : length(agglos)
    agglos_axon_reverse(agglos{idx}) = idx;
end
sizeThreshold = 1000:1000:300000;

for i=1:length(sizeThreshold)
    y(i).axon1 = connectEM.evaluateAggloMeta(skelpath(axon_selector), graph, segmentMeta, agglos_axon, p, 'axons1', mainFolder, 0, sizeThreshold(i), agglos_axon_reverse, 3000);
end
save('/gaba/scratch/mberning/coverage.mat', 'sizeThreshold', 'y');

