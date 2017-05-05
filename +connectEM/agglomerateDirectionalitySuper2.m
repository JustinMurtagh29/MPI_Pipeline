function agglomerateDirectionalitySuper(visualize)
 
    bboxDist = 1000;
    % Load needed data
    if visualize
        borderMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/globalBorder.mat', 'borderSize', 'borderCoM');
        graph = load('/gaba/u/mberning/results/pipeline/20170217_ROI/graphNew.mat', 'neighbours', 'neighProb', 'neighBorderIdx');
    else
        borderMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/globalBorder.mat', 'borderCoM');
        graph = load('/gaba/u/mberning/results/pipeline/20170217_ROI/graphNew.mat', 'neighbours', 'neighBorderIdx');
    end
    segmentMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta.mat', 'voxelCount', 'centroid', 'box', 'maxSegId');
    segmentMeta = connectEM.addSegmentClassInformation(p, segmentMeta);
    globalSegmentPCA = load('/gaba/u/mberning/results/pipeline/20170217_ROI/globalSegmentPCA.mat', 'covMat');
    gridAgglo_05{564} = load('/gaba/scratch/mberning/aggloGridSearch/search05_00564.mat', 'axonsFinal');
    % Agglomerates to use
    axonsFinalAll = cat(1, gridAgglo_05{564}.axonsFinal, ...
        num2cell(setdiff(find(segmentMeta.axonProb > 0.5), cell2mat(gridAgglo_05{564}.axonsFinal))));
    % Keep only agglomerates (including single segment agglomerados) over 100 voxel 
    minSize = 100;
    axonsFinalAllSize = cellfun(@(x)sum(segmentMeta.voxelCount(x)), axonsFinalAll);
    axonsFinalAll(axonsFinalAllSize < minSize) = [];
    y = connectEM.agglomerateDirectionality(axonsFinalAll(selection), graph, segmentMeta, borderMeta, globalSegmentPCA, bboxDist, visualize);

end

