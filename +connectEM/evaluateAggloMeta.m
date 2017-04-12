addpath(genpath('/gaba/u/kboerg/code/auxiliaryMethods'))
final = load('/gaba/scratch/mberning/20170404T172554_agglomerationNew/final.mat')
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameter.mat', 'p')
skel = skeleton('/gaba/u/kboerg/code/pipeline/+connectEM/evaluationData/2012-09-28_ex145_07x2_new2__explorational__mhelmstaedter__1bad47-2.nml')
segmentMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta.mat')
graph = load('/gaba/u/mberning/results/pipeline/20170217_ROI/graph.mat');
for idx = 1 : 2
    skel.nodes{idx}(:, 1 : 3) = bsxfun(@minus, skel.nodes{idx}(:, 1 : 3),  [1195, 1515, 115]-(129 - [25 25 10]));
    % skel = connectEM.evaluateAggloCleanSkel(skel, idx, any(skel.nodes{idx} <= 0, 2));
    ids = Seg.Global.getSegIds(p, skel.nodes{idx}(:, 1 : 3));
    skel = connectEM.evaluateAggloCleanSkel(skel, idx, ids == 0);
    skel = connectEM.intermediateNodes(skel, idx, 400)
    [splits, mergers] = connectEM.evaluateAgglo([final.final.axons; final.final.dendrites], segmentMeta.point', skel, idx, ids, graph.neighbours)
    splits_col(idx) = splits;
    mergers_col(idx) = mergers;
end
whos(matfile('/gaba/u/mberning/results/pipeline/20170217_ROI/graph.mat'))
