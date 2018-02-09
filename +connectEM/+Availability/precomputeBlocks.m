function [axonBlocks, targetClassAreas, targetClasses] = ...
        precomputeBlocks(param, conn, seg, blockSize)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    import connectEM.Availability.*;
    
    aggloSurfAreas = ...
        calculateSurfaceAreaOfAgglomerates(param, conn, seg, blockSize);
   [targetClassAreas, targetClasses] = ...
        calculateSurfaceAreaOfTargetClasses(conn, aggloSurfAreas);
    axonBlocks = calculateBlockListForAxons(conn, aggloSurfAreas);
end