function areas = axonSpineInterfaceArea( ...
        param, asiT, axons, spineHeads, varargin)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opts = struct;
    opts.box = [512; 512; 256];
    
    opts = Util.modifyStruct(opts, varargin{:});
    opts.box = reshape(opts.box, 3, 1);
    
    areas = nan(height(asiT), 4);
    for curId = 1:height(asiT)
        curAsi = asiT(curId, :);
        curAreas = cell(1, size(areas, 2));
        
       [curAreas{:}] = doIt( ...
            param, opts, asiT(curId, :), ...
            axons{curAsi.preAggloId}, spineHeads{curAsi.shId});
        curAreas = cell2mat(curAreas) %#ok
        areas(curId, :) = curAreas;
    end
end

function [areaRetina, areaSleep, areaIso, areaAlpha] = ...
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
        borders, edges, seg, param.raw.voxelSize, true);
    
    %% de Vivo et al. 2017 Science
    mask = false(size(seg));
    mask(borders.PixelIdxList) = true;
    
   [~, surfArea, capArea] = Seg.Local.physicalSurfaceArea( ...
        mask, param.raw.voxelSize, 'method', 'trakem2');
    areaSleep = (surfArea - capArea) / 2;
    areaSleep = areaSleep / (1E3) ^ 2;
    
    %% Isosurface
    iso = isosurface(mask, 0.5);
    iso.vertices = iso.vertices .* param.raw.voxelSize;
    iso = reducepatch(iso, 0.05);
    
    areaIso = cross( ...
        iso.vertices(iso.faces(:, 2), :) ...
      - iso.vertices(iso.faces(:, 1), :), ...
        iso.vertices(iso.faces(:, 3), :) ...
      - iso.vertices(iso.faces(:, 1), :), 2);
    areaIso = sum(sqrt(sum(areaIso .* areaIso, 2)) / 2);
    areaIso = areaIso / 2 / (1E3) ^ 2;
    
    %% Alpha Shape
    areaAlpha = ...
        Seg.Local.physicalBorderAreaAlphaShapes( ...
            borders, param.raw.voxelSize, size(seg));
    areaAlpha = areaAlpha / (1E3) ^ 2;
end
