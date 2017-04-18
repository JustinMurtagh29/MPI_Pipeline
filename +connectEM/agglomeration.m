function agglomeration( ... 
        borderSizeDendrites, segmentSizeDendrites, borderSizeAxons, segmentSizeAxons, ...
        axonProbThreshold, dendriteProbThreshold, spineProbThreshold, ...
        probThresholdDendrite, sizeThresholdDendrite, probThresholdAxon, sizeThresholdAxon, ...
        erProbThreshold, ...
        dendriteProbSpines, probThresholdSpines, maxStepsSpines, ...
        outputFile);

    %% Start by loading parameter file
    % Load parameter from newest pipeline run
    load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
    % To keep workspace clean here remove parameter for training (pT)
    clear pT;

    display('Loading data:');
    tic;
    % Load new graph representation including (resorted) borderIdx [p.saveFolder 'graphNew.mat']
    % These indicate which index in [p.saveFolder 'globalBorder.mat'] each edge corresponds to
    % Correspondences have NaN as borderIdx
    % Load 'neighbours' and 'neighProb' in addition if you want to do (many) local searches in the graph
    graph = load([p.saveFolder 'graphNew.mat'], 'prob', 'edges', 'borderIdx');
    % Load information about edges
    borderMeta = load([p.saveFolder 'globalBorder.mat'], 'borderSize', 'borderCoM');
    % Load meta information of segments
    segmentMeta = load([p.saveFolder 'segmentMeta.mat'], 'voxelCount', 'point', 'maxSegId');
    segmentMeta.point = segmentMeta.point';
    % Load and preprocess segment class predictions from Alessandro
    segmentMeta = connectEM.addSegmentClassInformation(p, segmentMeta);
    % Load exclusion of segments based on heuristics
    heuristics = load([p.saveFolder 'heuristicSegmentExclusions.mat'], 'mapping', 'heuristicIdx');
    toc;

    display('Removing segments detected by heuristics & small & disconnected segments:');
    tic;
    graphCutDendrites = connectEM.cutGraph(p, graph, segmentMeta, borderMeta, heuristics, ...
        borderSizeDendrites, segmentSizeDendrites);
    graphCutAxons = connectEM.cutGraph(p, graph, segmentMeta, borderMeta, heuristics, ...
        borderSizeAxons, segmentSizeAxons);
    toc;

    display('Generating subgraphs for axon and dendrite agglomeration:');
    tic;
    % Dendrites first
    idx = all(ismember(graphCutDendrites.edges, find(segmentMeta.dendriteProb > dendriteProbThreshold)), 2);
    graphCutDendrites.edges = graphCutDendrites.edges(idx,:);
    graphCutDendrites.prob = graphCutDendrites.prob(idx);
    % Then axon
    idx = all(ismember(graphCutAxons.edges, find(segmentMeta.axonProb > axonProbThreshold)), 2);
    graphCutAxons.edges = graphCutAxons.edges(idx,:);
    graphCutAxons.prob = graphCutAxons.prob(idx);
    clear idx;
    toc;

    display('Performing agglomeration on dendrite subgraph:');
    tic;
    [dendrites, dendriteSize, dendriteEdges] = connectEM.partitionSortAndKeepOnlyLarge(graphCutDendrites, segmentMeta, probThresholdDendrite, sizeThresholdDendrite);
    toc;

    display('Performing agglomeration on axon subgraph:');
    tic;
    [axons, axonsSize, axonEdges] = connectEM.partitionSortAndKeepOnlyLarge(graphCutAxons, segmentMeta, probThresholdAxon, sizeThresholdAxon);
    toc;
    
    %{
    display('Reassigning ER from axon to dendrite class: ');
    tic;
    [dendritesAfterEr, axonsAfterEr, er] = connectEM.extractAndTransferER(graph, dendrites, axons, erProbThreshold);
    toc;
    %}

    display('Garbage collection');
    tic;
    [axonsFinal, dendritesFinal] = connectEM.garbageCollection(graph, segmentMeta, axons, dendrites, heuristics.mapping);
    toc;

    %{
    display('Attaching spines to dendrite class: ');
    tic;
    [dendritesFinalWithSpines, spinePaths, comment] = connectEM.attachSpines(graph, segmentMeta, ...
        dendritesFinal, axonsFinal, spineProbThreshold, dendriteProbSpines, probThresholdSpines, maxStepsSpines);
    toc;

    display('Evaluating on a set of ground truth skeletons');
    tic;
    metrics = evalutateAggloMetaMeta(graph, cat(1, dendritesFinal, axonsFinal));  
    toc;
    %}

    display('Saving:');
    tic;
    % Lets save some more so that we can always understand whats happening
    save(outputFile, 'dendritesFinal', 'axonsFinal');
    toc;

end

