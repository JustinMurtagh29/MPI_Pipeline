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
    % Load global graph representation
    % graph = load([p.saveFolder 'graph.mat'], 'prob', 'edges', 'neighbours', 'neighProb'); 
    graph = load([p.saveFolder 'graph.mat'], 'prob', 'edges');
    % Load information about edges
    borderMeta = load([p.saveFolder 'globalBorder.mat'], 'borderSize', 'borderCoM');
    % Load meta information of segments
    segmentMeta = load([p.saveFolder 'segmentMeta.mat'], 'voxelCount', 'point', 'maxSegId');
    segmentMeta.point = segmentMeta.point';
    % Load and preprocess segment class predictions on single segments from Alessandro
    segmentMeta = connectEM.addSegmentClassInformation(p, segmentMeta);
    toc;

    display('Removing segments detected by heuristics & small & disconnected segments:');
    tic;
    [graphCutDendrites, excClassesDendrites, excNamesDendrites, excSizeDendrites] = connectEM.cutGraph(p, graph, segmentMeta, ...
        borderMeta, borderSizeDendrites, segmentSizeDendrites);
    graphCutAxons = connectEM.cutGraph(p, graph, segmentMeta, ...
        borderMeta, borderSizeAxons, segmentSizeAxons);
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
    % Note: Exlusions classes 1-6 consitent between dendrites and axons
    [axonsFinal, dendritesFinal] = connectEM.garbageCollection(graph, segmentMeta, axons, dendrites, excClassesDendrites(1:6));
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

