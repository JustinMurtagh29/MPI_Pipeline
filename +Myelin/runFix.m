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
    cluster = Cluster.config( ...
        'priority',500,...
        'memory',12, ...
        'time','12:00:00');
    job = Cluster.startJob(@Myelin.enforceMyelinSegments, taskInputArguments,'sharedInputs',{p}, 'cluster', cluster, 'name', 'myelinFix', 'taskGroupSize', 1);
end

