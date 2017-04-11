addpath(genpath('/gaba/u/kboerg/code/auxiliaryMethods'))
final = load('/gaba/scratch/mberning/20170404T172554_agglomerationNew/final.mat')
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameter.mat', 'p')
skel = skeleton('/gaba/u/kboerg/code/pipeline/+connectEM/evaluationData/2012-09-28_ex145_07x2_new2__explorational__mhelmstaedter__1bad47-2.nml')
segmentMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta.mat')
graph = load('/gaba/u/mberning/results/pipeline/20170217_ROI/graph.mat');
for idx = 1 : 2
    ids = Seg.Global.getSegIds(p, skel.nodes{idx}(:, 1 : 3))
    evaluateAgglo([final.axons; final.dendrites], segmentMeta.point, skel, idx, ids)
end
