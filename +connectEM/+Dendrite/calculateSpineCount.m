function spineCount = calculateSpineCount(param, dendrites, shAgglos)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    maxSegId = Seg.Global.getMaxSegId(param);
    shLUT = Agglo.buildLUT(maxSegId, shAgglos);
    
    spineCountFunc = @(segIds) numel(setdiff(shLUT(segIds), 0));
    spineCount = cellfun(spineCountFunc, dendrites);
end
