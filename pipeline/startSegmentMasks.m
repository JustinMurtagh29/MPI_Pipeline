function job = startSegmentMasks(param)
    cubes = param.local;
    cubeCount = numel(cubes);
    
    taskArgs = arrayfun( ...
        @(idx) {{cubes(idx)}}, 1:cubeCount);
    
    fH = @buildSegmentMasks;
    job = startCPU(fH, taskArgs, func2str(fH));
end