function job = aggloGridSearchDir()

    % axon related parameter in this run
    latentVar = [0.7 0.8 0.9]; % border below this size are excluded from axon subgraph during inital agglomeration
    segDirVar = [0.8 0.9 0.95]; % segment below ...
    neuriCVar = [0.6 0.7 0.8]; % restrictions of graph based on segment class proability
    borderVar = [20 40]; % threshold on neurite continuity probability for CC
    axonPrVar = [0.3 0.5]; % threshold on final agglomerate size in voxels (single segments not collected!)

    functionH = @connectEM.agglomerateDirectionalitySuper2;
    % Create parameter sets for axons and dendrites
    inputArgumentsAxons = inputArgumentsFromParameterSets( ...
        latentVar, segDirVar, neuriCVar, borderVar, axonPrVar);
    inputArgumentFolder = arrayfun(@(x)['/gaba/scratch/mberning/aggloGridSearch2/search01_' num2str(x, '%.5i') '/'], ...
        1:size(inputArgumentsAxons,1), 'uni', 0)';
    for i=1:size(inputArgumentsAxons,1)
        options.latentScore = inputArgumentsAxons(i,1);
        options.segDirScore = inputArgumentsAxons(i,2);
        options.neuriCScore = inputArgumentsAxons(i,3);
        options.borderSize = inputArgumentsAxons(i,4);
        options.axonScore = inputArgumentsAxons(i,5);
        inputArguments{i} = {options inputArgumentFolder{i}};
    end

    % Start job
    cluster = Cluster.getCluster( ... 
        '-pe openmp 1', ... 
        '-p -300', ...
        '-l h_vmem=48G', ... 
        '-l s_rt=99:50:00', ... 
        '-l h_rt=100:00:00');
    job = Cluster.startJob( functionH, inputArguments, ...
        'name', 'aggloGridSearchDir', ...
        'cluster', cluster);

end

function inputArguments = inputArgumentsFromParameterSets(a, b, c, d, e)

    [aG, bG, cG, dG, eG] = ndgrid(a, b, c, d, e);
    inputArguments = cat(2, aG(:), bG(:), cG(:), dG(:), eG(:)); 

end

