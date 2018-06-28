function [axonBlocks, targetClassAreas, targetClasses] = ...
        precompute(param, connFile, blockSize)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    %% per-cube calculations (on cluster)
    inputArgs = arrayfun(@(i) {{i}}, 1:numel(param.local));
    sharedInputArgs = {param, connFile, blockSize};

    job = Cluster.startJob( ...
        @precomputeCube, inputArgs, ...
        'sharedInputs', sharedInputArgs, ...
        'numOutputs', 3, ...
        'name', mfilename, ...
        'cluster', { ...
            'priority', 100, ...
            'memory', 24, ...
            'time', '1:00:00'});
    wait(job);
    
    %% assemblign output (locally)
    out = fetchOutputs(job);
    
    axonBlocks = cellfun( ...
        @vertcat, out{:, 1}, ...
        'UniformOutput', false);
    targetClasses = out{1, 3};
    
    blockCount = 1 + diff(param.bbox, 1, 2)';
    blockCount = ceil(blockCount ./ blockSize);
    
    targetClassAreas = nan(horzcat( ...
        numel(targetClasses), blockCount));
    for curIdx = 1:size(out, 1)
        curAreas = out{curIdx, 2};

        curOff = param.local(curIdx).bboxSmall(:, 1);
        curOff = (curOff - param.bbox(:, 1)) ./ blockSize(:);
        
        targetClassAreas(:, ...
            curOff(1) + (1:size(curAreas, 1 + 1)), ...
            curOff(2) + (1:size(curAreas, 1 + 2)), ...
            curOff(3) + (1:size(curAreas, 1 + 3))) = curAreas;
    end
end

function [axonBlocks, targetClassAreas, targetClasses] = ...
        precomputeCube(param, connFile, blockSize, cubeIdx)
    import connectEM.Availability.*;
    
    %% loading data
    conn = connectEM.Connectome.load(param, connFile);
    conn = connectEM.Connectome.prepareForSpecificityAnalysis(conn);
    
    box = param.local(cubeIdx).bboxSmall;
    seg = loadSegDataGlobal(param.seg, box);
    
    %% calculations
   [axonBlocks, targetClassAreas, targetClasses] = ...
       precomputeBlocks(param, conn, seg, blockSize);
   
    %% globalize block indices
    blockOff = param.local(cubeIdx).bboxSmall(:, 1);
    blockOff = (blockOff - param.bbox(:, 1)) ./ blockSize(:);
    blockOff = reshape(blockOff, 1, []);
    
    axonBlocks = cellfun( ...
        @(blocks) blocks + blockOff, ...
        axonBlocks, 'UniformOutput', false);
end
