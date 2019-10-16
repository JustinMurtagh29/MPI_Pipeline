function [borders, seg, box] = ...
        buildAxonSpineInterface(param, pos, boxSize, axonSegIds, shSegIds)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    invalidSegIds = intersect(axonSegIds, shSegIds);
    axonSegIds = setdiff(axonSegIds, invalidSegIds);
    shSegIds = setdiff(shSegIds, invalidSegIds);

    box = ceil(pos(:) + [-1, +1] .* boxSize(:) / 2);
    box = bsxfun(@max, box, param.bbox(:, 1));
    box = bsxfun(@min, box, param.bbox(:, 2));
    
    assert(all(box(:, 2) > box(:, 1)));
    seg = loadSegDataGlobal(param.seg, box);
    
    % For LUT
    seg = seg + 1;
    axonSegIds = axonSegIds + 1;
    shSegIds = shSegIds + 1;
    
    maxSegId = max(cellfun( ...
        @(ids) double(max(ids(:))), ...
        {seg, axonSegIds, shSegIds}));
    
    % NOTE(amotta): Map
    % * border voxels to label 0,
    % * axon segments to label 1,
    % * spine head segments to label 2, and
    % * all other segments to label 3.
    lut = zeros(maxSegId, 1, 'like', seg);
    lut(:) = 3; lut(1) = 0;
    lut(axonSegIds) = 1;
    lut(shSegIds) = 2;
    
    seg = lut(seg);
    seg(imclose(seg == 3, strel('cube', 3))) = 3;
    seg(imclose(seg == 2, strel('cube', 3))) = 2;
    seg(imclose(seg == 1, strel('cube', 3))) = 1;
    
    % NOTE(amotta): Closing doesn't work at end of cube.
    seg(1, :, :) = 3; seg(end, :, :) = 3;
    seg(:, 1, :) = 3; seg(:, end, :) = 3;
    seg(:, :, 1) = 3; seg(:, :, end) = 3;
    
   [edges, borders] = SynEM.Svg.findEdgesAndBorders(seg);
    borders = borders(all(edges == [1, 2], 2));
    assert(not(isempty(borders)));
end
