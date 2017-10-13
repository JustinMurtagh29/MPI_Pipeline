

graph=load('/gaba/u/mberning/results/pipeline/20170217_ROI/graphNewNew.mat')
borderMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/globalBorder.mat', 'borderSize', 'borderCoM');
segmentMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta.mat', 'voxelCount', 'point', 'maxSegId', 'cubeIdx', 'centroid', 'box');
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat', 'p');
segmentMeta = connectEM.addSegmentClassInformation(p, segmentMeta);
globalSegmentPCA = load('/gaba/u/mberning/results/pipeline/20170217_ROI/globalSegmentPCA.mat', 'covMat');
bboxDist = 5000;
[graph.neighbours, neighboursIdx] = Graph.edges2Neighbors(graph.edges);
graph.neighProb = cellfun(@(x)graph.prob(x), neighboursIdx, 'uni', 0);
graph.neighBorderIdx = cellfun(@(x)graph.borderIdx(x), neighboursIdx, 'uni', 0);
todo = [30    16    44    92    53    15    36     4    63    35    33    31    23     8    11    19    13    41    21    29    22    24    12     3    10    7];
parfor idx = 1:length(todo);
    idx
    y = connectEM.agglomerateDirectionality(dendritesFinal(todo(idx)), graph, segmentMeta, borderMeta, globalSegmentPCA, bboxDist, false)
    save(['/tmpscratch/kboerg/2017-10-12_directionalityDendrites' num2str(todo(idx)) '.mat'],'y')
end
