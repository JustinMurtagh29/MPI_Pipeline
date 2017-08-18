function lut = buildLUT(maxSegId, agglos)
    % lut = buildLUT(maxSegId, agglos)
    %   Builds a look-up table `lut` where entry `lut(segId)` contains the
    %   agglomerate ID the segment with ID `segId` belongs to. That is,
    %   
    %   ismember(segId, agglos{aggloId}) == true
    %   if and only if lut(segId) == aggloId
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    lut = zeros(maxSegId, 1);
    lut(cell2mat(agglos)) = ...
        cell2mat(arrayfun( ...
            @(idx) {repmat(idx, size(agglos{idx}))}, ...
            reshape(1:numel(agglos), size(agglos))));
end