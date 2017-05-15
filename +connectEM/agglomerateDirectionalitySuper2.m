function [axons, result] = agglomerateDirectionalitySuper2(options, outputFolder, graph, segmentMeta, borderMeta, globalSegmentPCA)

    % Create output folder if it does not exist
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder)
    end

    % Parameters
    recursionSteps = 10;
    minSize = 10;
    bboxDist = 1000;
    voxelSize = [11.24 11.24 28];

    % Load needed meta data if not passed (NOTE: this will take some time)
    if ~exist('options', 'var') | isempty(options)
        options.latentScore = 0.7; % 0.7-0.9
        options.segDirScore = 0.9; % 0.8-0.9
        options.neuriCScore = 0.7; % 0.6-0.9
        options.borderSize = 30; % 20:40
        options.axonScore = 0.3; % 0.3:0.5
    end
    if ~exist('graph', 'var') 
        graph = load('/gaba/u/mberning/results/pipeline/20170217_ROI/graphNew.mat', 'edges', 'prob', 'borderIdx');
        [graph.neighbours, neighboursIdx] = Graph.edges2Neighbors(edges);
        graph.neighProb = cellfun(@(x)prob(x), neighboursIdx, 'uni', 0);
        graph.neighBorderIdx = cellfun(@(x)borderIdx(x), neighboursIdx, 'uni', 0);
    end
    if ~exist('segmentMeta', 'var')
        segmentMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta.mat');
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
    [er, cm] = connectEM.getERcomponents();
    excludedCubeIdx = unique(cellfun(@(x)mode(segmentMeta.cubeIdx(x)), cm));
    excludedSegmentIdx = ismember(segmentMeta.cubeIdx, excludedCubeIdx) | ismember(1:segmentMeta.maxSegId, cat(1, er{:}))';
    % Add single segments if not excluded due to catastrohic merger
    axons = cat(1, gridAgglo_05{564}.axonsFinal, ...
        num2cell(setdiff(find(segmentMeta.axonProb > options.axonScore & ~excludedSegmentIdx), cell2mat(gridAgglo_05{564}.axonsFinal))));
    clear gridAgglo_05 excluded* cm;

    % Keep only agglomerates (including single segment agglomerados) over minSize voxel 
    axonsSize = cellfun(@(x)sum(segmentMeta.voxelCount(x)), axons);
    axons(axonsSize < minSize) = [];
    % For debugging: Keep only 10000 random components
    axons = axons(randperm(numel(axons), 10000));

    % Initialize variables that keep track of changed agglomerates over recursion for first round 
    changedIdx = 1:length(axons);
    unchangedResult = struct('latent', [], 'pca', [], 'neighbours', [], 'prob', [], 'borderIdx', [], 'scores', []);
    axonsNew = axons;

    for i=1:recursionSteps
        axons = axonsNew;
        clear axonsNew result;
        display('Calculating directionality measures:');
        tic;
        resultTemp = connectEM.agglomerateDirectionality2(axons(changedIdx), graph, segmentMeta, borderMeta, globalSegmentPCA, bboxDist, voxelSize);
        fieldNames = fieldnames(resultTemp);
        for j=1:length(fieldNames)
            fName = fieldNames{j};
            result.(fName) = cat(1, resultTemp.(fName), unchangedResult.(fName));
        end
        clear resultTemp fieldNames j changedIdx unchangedResult;
        toc;
        display('Merging agglomerates:');
        tic;
        [axonsNew, changedIdx, unchangedResult] = connectEM.agglomerateMerge(graph, segmentMeta, borderMeta, axons, result, options);
        toc;
        display('Saving (intermediate) results:');
        tic;
        Util.save([outputFolder num2str(i, '%.2i') '.mat'], result, axons, axonsNew);
        toc;
    end

end

