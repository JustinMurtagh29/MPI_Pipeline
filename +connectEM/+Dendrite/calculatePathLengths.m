function pathLens = calculatePathLengths(param, dendrites, trunkAgglos)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    maxSegId = Seg.Global.getMaxSegId(param);
    segPoints = Seg.Global.getSegToPointMap(param);
    
    %% Reduce to trunk
    mask = false(maxSegId, 1);
    mask(cell2mat(trunkAgglos)) = true;

    dendrites = cellfun( ...
        @(segIds) segIds(mask(segIds)), ...
        dendrites, 'UniformOutput', false);

    dendrites = cellfun( ...
        @(segIds) segPoints(segIds, :), ...
        dendrites, 'UniformOutput', false);
    dendrites = struct('nodes', dendrites);
    
    %% Calculate path length
    pathLens = Superagglos.mstLength( ...
        dendrites, param.raw.voxelSize);
end
