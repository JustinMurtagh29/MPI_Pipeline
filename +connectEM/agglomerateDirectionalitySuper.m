function agglomerateDirectionalitySuper(slice)

visualize = false;
if visualize
    graph = load('/gaba/u/mberning/results/pipeline/20170217_ROI/graphNew.mat', 'neighbours', 'neighProb', 'neighBordIdx')
else
    graph = load('/gaba/u/mberning/results/pipeline/20170217_ROI/graphNew.mat', 'neighbours', 'neighBordIdx')
end
segmentMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta.mat', 'voxelCount', 'centroid', 'box');
globalSegmentPCA = load('/gaba/u/mberning/results/pipeline/20170217_ROI/globalSegmentPCA.mat', 'covMat');
borderMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/borderMeta.mat');
gridAgglo_05{564} = load('/gaba/scratch/mberning/aggloGridSearch/search05_00564.mat', 'axonsFinal');

bboxDist = 1000;
axonsFinalAll = [gridAgglo_05{564}.axonsFinal,
    num2cell(setdiff(1 : length(segmentMeta.voxelCount), cell2mat(gridAgglo_05{564}.axonsFinal)))'];
selection = (1 + 50000 * (slice - 1)) : min(50000 * slice, length(axonsFinalAll));
y = connectEM.agglomerateDirectionality(axonsFinalAll(selection), graph, segmentMeta, borderMeta, globalSegmentPCA, bboxDist, visualize)
save(['/gaba/scratch/kboerg/cluster_directionality_run/', num2str(slice, '%.4u') '.mat'], 'y')
