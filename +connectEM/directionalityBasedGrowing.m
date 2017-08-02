function directionalityBasedGrowing(options, outputFolder, agglos, graph, segmentMeta, borderMeta, globalSegmentPCA, heuristics)

    % Create output folder if it does not exist
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder)
    end

    if ~exist('options', 'var') || isempty(options)
        options.latentScore = 0.7;
        options.segDirScore = 0.9;
        options.neuriCScore = 0.7;
        options.borderSize = 30;
        options.axonScore = 0.3;
        options.sourceSize = 2000;
        options.recursionSteps = 10;
        options.minSize = 100;
        options.bboxDist = 1000;
        options.voxelSize = [11.24 11.24 28];
    end

    % Load needed meta data if not passed
    if ~exist('graph', 'var')
        graph = load('/gaba/u/mberning/results/pipeline/20170217_ROI/graphNew.mat', 'edges', 'prob', 'borderIdx');
        [graph.neighbours, neighboursIdx] = Graph.edges2Neighbors(graph.edges);
        graph.neighProb = cellfun(@(x)graph.prob(x), neighboursIdx, 'uni', 0);
        graph.neighBorderIdx = cellfun(@(x)graph.borderIdx(x), neighboursIdx, 'uni', 0);
    end
    if ~exist('segmentMeta', 'var')
        load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
        segmentMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta.mat');
        segmentMeta = connectEM.addSegmentClassInformation(p, segmentMeta);
    end
    if ~exist('borderMeta', 'var')
        borderMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/globalBorder.mat', 'borderSize', 'borderCoM');
    end
    if ~exist('globalSegmentPCA', 'var')
        globalSegmentPCA = load('/gaba/u/mberning/results/pipeline/20170217_ROI/globalSegmentPCA.mat', 'covMat');
    end
    if ~exist('heuristics', 'var')
        heuristics = load('/gaba/u/mberning/results/pipeline/20170217_ROI/heuristicResult.mat');
    end
%     % Exclude all segments in cubes with catastrophic merger
%     [er, cm] = connectEM.getERcomponents();
%     excludedCubeIdx = unique(cellfun(@(x)mode(segmentMeta.cubeIdx(x)), cm));
%     excludedSegmentIdx = ismember(segmentMeta.cubeIdx, excludedCubeIdx) | ismember(1:segmentMeta.maxSegId, cat(1, er{:}))';
%     % Add single segments if not excluded due to catastrophic merger or ER 
% %     startAgglo = load('/gaba/scratch/mberning/aggloGridSearch/search05_00564.mat', 'axonsFinal');
% %     agglos = cat(1, startAgglo.axonsFinal, ...
% %         num2cell(setdiff(find((segmentMeta.axonProb > options.axonScore | segmentMeta.dendriteProb > options.dendriteScore) & ~excludedSegmentIdx), cell2mat(startAgglo.axonsFinal))));
% %     agglos = num2cell(find((segmentMeta.axonProb > options.axonScore | segmentMeta.dendriteProb > options.dendriteScore) & ~excludedSegmentIdx));
%     clear startAgglo excluded* cm er;

    % Keep only agglomerates (including single segment agglomerados) over minSize voxel
    aggloSize = cellfun(@(x)sum(segmentMeta.voxelCount(x)), agglos);
    agglos(aggloSize < options.minSize) = [];

    % Initialize variables that keep track of changed agglomerates over recursion for first round
    changedIdx = 1:length(agglos);
    unchangedResult = struct('latent', [], 'pca', [], 'neighbours', [], 'prob', [], 'borderIdx', [], 'scores', []);
    agglosNew = agglos;

    for i=1:options.recursionSteps
        
        agglos = agglosNew;
        clear axonsNew result;
        
        display('Calculating directionality measures:');
        tic;
        resultTemp = connectEM.calculateDirectionalityOfAgglomerates(agglos(changedIdx), graph, segmentMeta, borderMeta, globalSegmentPCA, options);
        fieldNames = fieldnames(resultTemp);
        for j=1:length(fieldNames)
            fName = fieldNames{j};
            result.(fName) = cat(1, resultTemp.(fName), unchangedResult.(fName));
        end
        clear resultTemp fieldNames j changedIdx unchangedResult;
        toc;
        
        myelinScore = cellfun(@(x) heuristics.myelinScore(x) < 0.5,agglos);  % take only segments with a low myelinScore

        display('Merging agglomerates:');
        tic;
        [agglosNew, changedIdx, unchangedResult, edgesToStore] = connectEM.agglomerateMerge(graph, segmentMeta, borderMeta, agglos, result, options, myelinScore);
        toc;
        if options.recursionSteps > 1
            display('Saving (intermediate) results:');
            tic;
            Util.save([outputFolder num2str(i, '%.2i') '.mat'], agglosNew, edgesToStore);
            toc;
        end
    end
    
    display('Calculating set of extended metrics:');
    tic;
    name = strsplit(outputFolder, '/');
    name = name{end-1};
    metrics = connectEM.moreMetrics(agglosNew, name, segmentMeta);
    toc;

    display('Saving metrics:');
    tic;
    Util.save([outputFolder 'metricsFinal.mat'], agglosNew, metrics);
    toc;

end

