function areas = buildAxonSpineInterfaceAreas( ...
        param, axons, spineHeads, asiT, varargin)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opts = struct;
    opts.box = [512; 512; 256];
    
    opts = Util.modifyStruct(opts, varargin{:});
    opts.box = reshape(opts.box, 3, 1);
    
    areas = nan(height(asiT), 1);
    for curId = 1:height(asiT)
        curAsi = asiT(curId, :);
        
        areas(curId, :) = doIt( ...
            param, opts, curAsi.pos, ...
            axons{curAsi.preAggloId}, ...
            spineHeads{curAsi.shId});
    end
end

function area = doIt(param, opts, pos, axonSegIds, shSegIds)
   [borders, seg] = ...
        connectEM.Consistency.buildAxonSpineInterfaces( ...
            param, pos, opts.box, axonSegIds, shSegIds);
    area = ...
        Seg.Local.physicalBorderAreaAlphaShapes( ...
            borders, param.raw.voxelSize, size(seg));
    area = sum(area) / 1E6;
end
