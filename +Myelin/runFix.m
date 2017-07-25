function job = runFix(p, newPrefix)
    % job = runFix(p, newPrefix)
    %   Applies `+Myelin/runFixInBox` to the entire dataset.
    %   For further information, check out the documentation
    %   of said function.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    boxCount = numel(p.local);
    taskInputArguments = arrayfun( ...
        @(idx) {{newPrefix, p.local(idx).bboxSmall}}, 1:boxCount);
    %     job = startCPU(@Myelin.enforceMyelinSegments, taskInputArguments, 'myelinFix');
    cluster = Cluster.getCluster( ...
        '-pe openmp 1', ...
        '-p -500', ...
        '-l h_vmem=12G', ...
        '-l s_rt=12:00:00', ...
        '-l h_rt=12:00:30');
    job = Cluster.startJob(@Myelin.enforceMyelinSegments, taskInputArguments,'sharedInputs',{p}, 'cluster', cluster, 'name', 'myelinFix', 'taskGroupSize', 1);
end

