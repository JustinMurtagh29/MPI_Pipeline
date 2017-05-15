function [axons, result] = agglomerateDirectionalitySuper2(graph, segmentMeta, borderMeta, globalSegmentPCA)
 
    % Parameters
    recursionSteps = 3;
    minSize = 100;
    bboxDist = 1000;
    voxelSize = [11.24 11.24 28];
    
    % Load needed meta data if not passed (NOTE: this will take some time)
    if ~exist('graph', 'var') 
        graph = load('/gaba/u/mberning/results/pipeline/20170217_ROI/graphNew.mat');
    end
    if ~exist('segmentMeta', 'var')
        segmentMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta.mat', 'voxelCount', 'centroid', 'box', 'maxSegId', 'cubeIdx');
        segmentMeta = connectEM.addSegmentClassInformation(p, segmentMeta);
    end
    if ~exist('borderMeta', 'var')
        borderMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/globalBorder.mat', 'borderSize', 'borderCoM');
    end 
    if ~exist('globalSegmentPCA', 'var')
        globalSegmentPCA = load('/gaba/u/mberning/results/pipeline/20170217_ROI/globalSegmentPCA.mat', 'covMat');
    end
    
    % Agglomerates to start with
    gridAgglo_05{564} = load('/gaba/scratch/mberning/aggloGridSearch/search05_00564.mat', 'axonsFinal');
    % Exclude all segments in cubes with catastrophic merger
    [~, cm] = connectEM.getERcomponents();
    excludedCubeIdx = unique(cellfun(@(x)mode(segmentMeta.cubeIdx(x)), cm));
    excludedSegmentIdx = ismember(segmentMeta.cubeIdx, excludedCubeIdx);
    % Add single segments if not excluded due to catastrohic merger
    axons{1} = cat(1, gridAgglo_05{564}.axonsFinal, ...
        num2cell(setdiff(find(segmentMeta.axonProb > 0.3 & ~excludedSegmentIdx), cell2mat(gridAgglo_05{564}.axonsFinal))));
    clear gridAgglo_05 excluded* cm;
    
    % Keep only agglomerates (including single segment agglomerados) over minSize voxel 
    axonsSize = cellfun(@(x)sum(segmentMeta.voxelCount(x)), axons{1});
    axons{1}(axonsSize < minSize) = [];
    % For debugging: Keep only 10000 random components
    %axons{1} = axons{1}(randperm(numel(axons{1}), 10000));

    for i=1:recursionSteps
        display('Calculating directionality measures:');
        tic;
        result(i) = connectEM.agglomerateDirectionality2(axons{i}, graph, segmentMeta, borderMeta, globalSegmentPCA, bboxDist, voxelSize);
        toc;
        display('Merging agglomerates:');
        tic;
        axons{i+1} = connectEM.agglomerateMerge(graph, segmentMeta, borderMeta, axons{i}, result(i));
        toc;
        display('Saving (intermediate) results:');
        tic;
        Utilsave(['/gaba/scratch/mberning/dirTest/' num2str(i, '%.2i') '.mat'], result, axons, '-v7.3');
        toc;
    end

end

