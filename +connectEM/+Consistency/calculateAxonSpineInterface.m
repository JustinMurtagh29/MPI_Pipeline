function areas = calculateAxonSpineInterface( ...
        param, asiT, axons, spineHeads, varargin)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opts = struct;
    opts.box = [512; 512; 256];
    
    opts = Util.modifyStruct(opts, varargin{:});
    opts.box = reshape(opts.box, 3, 1);
    
    areas = nan(height(asiT), 2);
    for curId = 1:height(asiT)
        curAsi = asiT(curId, :);

       [areas(curId, 1), areas(curId, 2)] = doIt( ...
            param, opts, asiT(curId, :), ...
            axons{curAsi.preAggloId}, spineHeads{curAsi.shId});
    end
end

function [areaRetina, areaSleep] = ...
        doIt(param, opts, asi, axonSegIds, shSegIds)
    invalidSegIds = intersect(axonSegIds, shSegIds);
    axonSegIds = setdiff(axonSegIds, invalidSegIds);
    shSegIds = setdiff(shSegIds, invalidSegIds);

    box = ceil(asi.pos(:) + [-1, +1] / 2 .* opts.box);
    seg = loadSegDataGlobal(param.seg, box);
    
    % For LUT
    seg = seg + 1;
    axonSegIds = axonSegIds + 1;
    shSegIds = shSegIds + 1;
    
    maxSegId = max(cellfun( ...
        @(ids) double(max(ids(:))), ...
        {seg, axonSegIds, shSegIds}));
    
    lut = zeros(maxSegId, 1, 'like', seg);
    lut(:) = 3; lut(1) = 0;
    lut(axonSegIds) = 1;
    lut(shSegIds) = 2;
    
    seg = lut(seg);
    seg(imclose(seg == 3, strel('cube', 3))) = 3;
    seg(imclose(seg == 2, strel('cube', 3))) = 2;
    seg(imclose(seg == 1, strel('cube', 3))) = 1;
    
    % NOTE(amotta): Closing doesn't work at end of cube.
    seg(1, :, :) = 0; seg(end, :, :) = 0;
    seg(:, 1, :) = 0; seg(:, end, :) = 0;
    seg(:, :, 1) = 0; seg(:, :, end) = 0;
    
    %% Helmstaedter et al. 2013 Nature
   [edges, borders] = ...
        SynEM.Svg.findEdgesAndBorders(seg);
    
   [~, mask] = ismember([1, 2], edges, 'rows');
    edges = edges(mask, :);
    borders = borders(mask);
    
    areaRetina = Seg.Local.physicalBorderArea2( ...
        borders, edges, seg, param.raw.voxelSize, true) %#ok
    
    %% de Vivo et al. 2017 Science
    % via surface area measurement method of TrakEM2 
    mask = false(size(seg));
    mask(borders.PixelIdxList) = true;
    
   [~, curSurfArea, curCapArea] = ...
        Seg.Local.physicalSurfaceArea(mask, param.raw.voxelSize);
    areaSleep = (curSurfArea - curCapArea) / 2 / (1E3) ^ 2 %#ok
end
