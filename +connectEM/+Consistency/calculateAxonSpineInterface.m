function areas = calculateAxonSpineInterface( ...
        param, asiT, axons, spineHeads, varargin)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opts = struct;
    opts.box = [512; 512; 256];
    
    opts = Util.modifyStruct(opts, varargin{:});
    opts.box = reshape(opts.box, 3, 1);
    
    areas = nan(height(asiT), 1);
    for curId = 1:height(asiT)
        curAsi = asiT(curId, :);

        areas(curId) = doIt( ...
            param, opts, asiT(curId, :), ...
            axons{curAsi.preAggloId}, spineHeads{curAsi.shId});
    end
end

function areaSleep = doIt(param, opts, asi, axonSegIds, shSegIds)
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
    mask = false(size(seg));
    mask(borders.PixelIdxList) = true;
    
    maskIdsZ = mask;
    maskIdsZ = shiftdim(any(maskIdsZ, 1), 1);
    maskIdsZ = shiftdim(any(maskIdsZ, 1), 1);
    maskIdsZ = reshape(find(maskIdsZ), 1, []);
    
    areas = zeros(size(mask, 3), 1);
    perim = zeros(size(mask, 3), 1);
    
    for curId = maskIdsZ
        curMask = mask(:, :, curId);
        areas(curId) = bwarea(curMask);
        
        curPerim = regionprops(curMask, 'Perimeter');
        perim(curId) = curPerim.Perimeter;
    end
    
    areas = areas .* prod(param.raw.voxelSize(1:2));
    perim = perim .* param.raw.voxelSize(1);
    
    areaSleep = ...
        perim(1:(end - 1)) * param.raw.voxelSize(3) / 2 ...
      + perim(2:end) * param.raw.voxelSize(3) / 2 ...
      + areas(1:(end - 1)) - areas(2:end);
    areaSleep = sum(areaSleep);
  
    % Correct for capping at top and bottom
    areaSleep = ...
        areaSleep ...
      - areas(maskIdsZ(1)) ...
      - areas(maskIdsZ(end));
  
    % Correct for two-facedness
    areaSleep = areaSleep / 2;
    areaSleep = areaSleep / (1E3 ^ 2) %#ok
end
