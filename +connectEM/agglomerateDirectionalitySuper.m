function agglomerateDirectionalitySuper(slice, visualize, startagglo, axonProbThreshold, outputfolder)

if visualize
    borderMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/globalBorder.mat', 'borderSize', 'borderCoM');
    graph = load('/gaba/u/mberning/results/pipeline/20170217_ROI/graphNew.mat', 'neighbours', 'neighProb', 'neighBorderIdx');
else
    borderMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/globalBorder.mat', 'borderCoM');
    graph = load('/gaba/u/mberning/results/pipeline/20170217_ROI/graphNew.mat', 'neighbours', 'neighBorderIdx');
end
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat', 'p');
segmentMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta.mat', 'voxelCount', 'point', 'maxSegId', 'cubeIdx', 'centroid', 'box');
segmentMeta = connectEM.addSegmentClassInformation(p, segmentMeta);

globalSegmentPCA = load('/gaba/u/mberning/results/pipeline/20170217_ROI/globalSegmentPCA.mat', 'covMat');
agglo_col = load(startagglo, 'axonsFinal');
bboxDist = 1000;
axonsFinalAll = [agglo_col.axonsFinal,
    num2cell(setdiff(find(segmentMeta.axonProb > axonProbThreshold), cell2mat(agglo_col.axonsFinal)))];
selection = slice:200:length(axonsFinalAll);
y = connectEM.agglomerateDirectionality(axonsFinalAll(selection), graph, segmentMeta, borderMeta, globalSegmentPCA, bboxDist, visualize);
Util.save([outputfolder, num2str(slice, '%.4u') '.mat'], y);

end
