function removeOverlaps(p)
    % Transfer local cubes to segFile
    inputCell = arrayfun( ...
        @(local) { ...
            local.tempSegFile, local.segFile, ...
            local.bboxSmall, local.bboxBig}, ...
        p.local(:), 'UniformOutput', false);
    
    functionH = @rewriteWithoutOverlaps;
    job = startCPU(functionH, inputCell, 'overlapRemoval', 12, 10);
    Cluster.waitForJob(job);
    rmdir(fileparts(fileparts(p.local(1).tempSegFile)),'s') % remove temp structure (which is huge)
end

