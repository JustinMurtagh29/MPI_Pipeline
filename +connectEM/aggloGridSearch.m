function job = aggloGridSearch()

    % dendrite related parameter in this run
    borderSizeDendrites = [50 100 200 300 500]; % border below this size are excluded from dendrite subgraph during inital agglomeration
    segmentSizeDendrites = [100 300 500 700 1000]; % segment below ...
    dendriteProbThreshold = [0 0.1 0.2 0.3 0.4 0.5]; % restrictions of graph based on segment class proability
    probThresholdDendrite = [0.94 0.95 0.96 0.97 0.98 0.99]; % threshold on neurite continuity probability for CC
    sizeThresholdDendrite = 100; % threshold on final agglomerate size in voxels (single segments not collected!)

    % axon related parameter in this run
    borderSizeAxons = [0 25 50 100 200]; % border below this size are excluded from axon subgraph during inital agglomeration
    segmentSizeAxons = [50 100 200 300 400]; % segment below ...
    axonProbThreshold = [0 0.1 0.2 0.3 0.4 0.5]; % restrictions of graph based on segment class proability
    probThresholdAxon = [0.90 0.91 0.92 0.93 0.94 0.95]; % threshold on neurite continuity probability for CC
    sizeThresholdAxon = 100; % threshold on final agglomerate size in voxels (single segments not collected!)

    % NOTE: Parameters currently have no effect, commented routines for first tests
    % ... for ER reassignment
    erProbThreshold = 2;
    % ... for spine attachment
    spineProbThreshold = 0.5;
    dendriteProbSpines = 0.05;
    probThresholdSpines = 0.05;
    maxStepsSpines = 10;

    functionH = @connectEM.agglomeration;
    % Create parameter sets for axons and dendrites
    inputArgumentsDendrites = inputArgumentsFromParameterSets( ...
        borderSizeDendrites, segmentSizeDendrites, dendriteProbThreshold, probThresholdDendrite, sizeThresholdDendrite);
    inputArgumentsAxons = inputArgumentsFromParameterSets( ...
        borderSizeAxons, segmentSizeAxons, axonProbThreshold, probThresholdAxon, sizeThresholdAxon);
    % Currently 'static' parameter and filenames to save results to
    inputArgumentsStatic = repmat([erProbThreshold spineProbThreshold dendriteProbSpines probThresholdSpines maxStepsSpines], ...
        size(inputArgumentsDendrites, 1), 1);
    inputArgumentFilename = arrayfun(@(x)['/gaba/scratch/mberning/aggloGridSearch/search01_' num2str(x, '%.5i') '.mat'], ...
        1:size(inputArgumentsDendrites,1), 'uni', 0)';
    % Collect and reorder input arguments together
    inputArguments = cat(2, inputArgumentsDendrites(:,1:2), inputArgumentsAxons(:,1:2), ...
        inputArgumentsAxons(:,3), inputArgumentsDendrites(:,3), inputArgumentsStatic(:,2), ...
        inputArgumentsDendrites(:,4:5), inputArgumentsAxons(:,4:5), inputArgumentsStatic(:,1), ...
        inputArgumentsStatic(:,3:5));
    inputArguments = mat2cell(inputArguments, ones(size(inputArguments,1),1), size(inputArguments,2));
    inputArguments = cellfun(@(x,y)[num2cell(x) y], inputArguments, inputArgumentFilename, 'uni', 0);

    % Start job
    cluster = Cluster.getCluster( ... 
                '-pe openmp 1', ... 
                '-p 0', ...
                '-l h_vmem=12G', ... 
                '-l s_rt=23:50:00', ... 
                '-l h_rt=24:00:00');
    job = Cluster.startJob( functionH, inputArguments, ...
                'name', 'aggloGridSearch', ...
                'cluster', cluster);

end

function inputArguments = inputArgumentsFromParameterSets(a, b, c, d, e)

    [aG, bG, cG, dG, eG] = ndgrid(a, b, c, d, e);
    inputArguments = cat(2, aG(:), bG(:), cG(:), dG(:), eG(:)); 

end

