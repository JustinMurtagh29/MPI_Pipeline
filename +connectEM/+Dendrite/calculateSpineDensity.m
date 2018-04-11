function [spineDensity, spineCount] = ...
        calculateSpineDensity(param, dendrites, dendLengths, shAgglos)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    maxSegId = Seg.Global.getMaxSegId(param);
    shLUT = Agglo.buildLUT(maxSegId, shAgglos);
    
    spineIds = cellfun( ...
        @(segIds) setdiff(shLUT(segIds), 0), ...
        dendrites, 'UniformOutput', false);
    
    spineCount = cellfun(@numel, spineIds);
    spineDensity = spineCount ./ dendLengths;
end
