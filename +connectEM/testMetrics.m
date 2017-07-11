
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
evalfolder = '/gaba/u/kboerg/code/pipeline/+connectEM/evaluationData/';
% Load one axon GT and lookup segIds as test case
skel = skeleton([evalfolder '2012-09-28_ex145_07x2_new2__explorational__mhelmstaedter__1bad47-2_axon.nml']);
skel.nodes{1}(:, 1 : 3) = bsxfun(@minus, skel.nodes{1}(:, 1 : 3),  [1195, 1515, 115]-(129 - [25 25 10]));
skel.nodes{1}(skel.nodes{1} <= 0) = 1;
skel.nodesNumDataAll{1}(:, 3 : 5)= skel.nodes{1}(:, 1 : 3);
segIdsA = Seg.Global.getSegIds(p, skel.nodes{1}(:, 1 : 3));
segIdsA(segIdsA == 0) = [];
% Load one dendrite GT and lookup segIds as test case
skel = skeleton([evalfolder '2012-09-28_ex145_07x2_new2__explorational__kboergens__bd3e11.nml']);
skel.nodes{1}(:, 1 : 3) = bsxfun(@minus, skel.nodes{1}(:, 1 : 3),  [1195, 1515, 115]-(129 - [25 25 10]));
skel.nodes{1}(skel.nodes{1} <= 0) = 1;
skel.nodesNumDataAll{1}(:, 3 : 5)= skel.nodes{1}(:, 1 : 3);
segIdsD = Seg.Global.getSegIds(p, skel.nodes{1}(:, 1 : 3));
segIdsD(segIdsD == 0) = [];

segmentMeta = load([p.saveFolder 'segmentMeta.mat'], 'voxelCount', 'point', 'maxSegId');
segmentMeta.point = segmentMeta.point';
y = connectEM.evaluateAggloMetaMeta(struct('neighbours', 0), {segIdsA}, {segIdsD}, 'test345', segmentMeta);

