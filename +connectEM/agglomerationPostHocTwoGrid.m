function job = aggloPostHocTwoGrid()

    % axon related parameter in this run
    optionsPre.borderSizeThreshold = [30 40 50 60 70]; % border below this size are excluded from axon subgraph during inital agglomeration
    optionsPre.sizeThreshold = [0 100 200 250 300 400 500]; % segment below ...
    optionsPre.axonProbThreshold = [0.1 0.2 0.3 0.4 0.5]; % restrictions of graph based on segment class proability
    optionsPre.probThreshold = [0.93 0.94 0.95 0.96 0.97 0.98]; % threshold on neurite continuity probability for CC
    optionsPre.scoreThreshold = [0.8 0.85 0.9 0.95 0.98];
    optionsPre.latentThreshold = [0.6 0.7 0.8];
    optionsPre.agglomerationSizeThreshold = [100, 1000, 2500, 6250]; % threshold on final agglomerate size in voxels (single segments not collected!)


    functionH = @connectEM.agglomerationPostHocTwo;
    inputArgumentsAxons = inputArgumentsFromParameterSets(optionsPre);
    inputArgumentFilename = arrayfun(@(x)['/gaba/scratch/mberning/aggloGridSearch/searchPostHoc01_' num2str(x, '%.5i') '.mat'], ...
        1:size(inputArgumentsAxons,1), 'uni', 0)';
    % Collect and reorder input arguments together
    inputArguments = cellfun(@(x,y)[x y], inputArguments, inputArgumentFilename, 'uni', 0);

    % Start job
    cluster = Cluster.getCluster( ...
                '-pe openmp 1', ...
                '-p -400', ...
                '-l h_vmem=23G', ...
                '-l s_rt=23:50:00', ...
                '-l h_rt=24:00:00');
    job = Cluster.startJob( functionH, inputArguments, ...
                'name', 'aggloGridSearch', ...
                'cluster', cluster);

end

function inputArguments = inputArgumentsFromParameterSets(optionsPre)
    fns = fieldnames(optionsPre);
    for idx = 1 : length(fns)
        todo
    [aG, bG, cG, dG, eG] = ndgrid(a, b, c, d, e);
    inputArguments = cat(2, aG(:), bG(:), cG(:), dG(:), eG(:));

end
