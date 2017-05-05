function agglomerateDirectionalitySuper(slice, visualize)

if visualize
    borderMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/globalBorder.mat', 'borderSize', 'borderCoM');
    graph = load('/gaba/u/mberning/results/pipeline/20170217_ROI/graphNew.mat', 'neighbours', 'neighProb', 'neighBorderIdx');
else
    borderMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/globalBorder.mat', 'borderCoM');
    graph = load('/gaba/u/mberning/results/pipeline/20170217_ROI/graphNew.mat', 'neighbours', 'neighBorderIdx');
end
segmentMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta.mat', 'voxelCount', 'centroid', 'box');
globalSegmentPCA = load('/gaba/u/mberning/results/pipeline/20170217_ROI/globalSegmentPCA.mat', 'covMat');
gridAgglo_05{564} = load('/gaba/scratch/mberning/aggloGridSearch/search05_00564.mat', 'axonsFinal');
bboxDist = 1000;
axonsFinalAll = [gridAgglo_05{564}.axonsFinal,
    num2cell(setdiff(1 : length(segmentMeta.voxelCount), cell2mat(gridAgglo_05{564}.axonsFinal)))'];
selection = slice:200:length(axonsFinalAll);
y = connectEM.agglomerateDirectionality(axonsFinalAll(selection), graph, segmentMeta, borderMeta, globalSegmentPCA, bboxDist, visualize);
Util.save(['/gaba/scratch/kboerg/cluster_directionality_run2/', num2str(slice, '%.4u') '.mat'], y);

end
