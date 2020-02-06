function buildSynapseIsosurfaces( ...
        param, axons, dendrites, outDir, synT, varargin)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    info = Util.runInfo(false);
    
    opts = struct;
    opts.box = [512; 512; 256];
    opts.sphereRadNm = 50;
    
    opts = Util.modifyStruct(opts, varargin{:});
    opts.box = reshape(opts.box, 3, 1);
    opts.info = info;
    
    for curId = 1:height(synT)
        curSyn = synT(curId, :);
        
        curOutFile = sprintf('iso-%d.mat', curSyn.id);
        curOutFile = fullfile(outDir, curOutFile);
        
        doIt( ...
            param, opts, curSyn.pos, ...
            axons{curSyn.preAggloId}, ...
            dendrites{curSyn.postAggloId}, ...
            curOutFile);
    end
end

function isoSurf = doIt( ...
        param, opts, pos, axonSegIds, dendSegIds, outFile)
   [borders, seg, box] = ...
        connectEM.Consistency.buildAxonSpineInterface( ...
            param, pos, opts.box, axonSegIds, dendSegIds);
    borders = cell2mat(arrayfun( ...
        @(b) b.PixelIdxList, borders(:), ...
        'UniformOutput', false));
    
    sphereRadNm = opts.sphereRadNm;
   [~, sphereRadVx] = Util.ballIdx( ...
       sphereRadNm, size(seg), param.raw.voxelSize);
    
    assert(all(mod(sphereRadVx, 2) == 1));
    sphereRadVx = (sphereRadVx - 1) / 2;
    
    sizeRaw = size(seg);
    sizePadded = sizeRaw + 2 * sphereRadVx;
    
   [sphereOffsets, ~] = Util.ballIdx( ...
       sphereRadNm, sizePadded, param.raw.voxelSize);
    
    borders = Util.indToSubMat(sizeRaw, borders);
    borders = Util.sub2ind(sizePadded, borders + sphereRadVx);
    borders = borders(:) + transpose(sphereOffsets(:));
    borders = unique(borders(:));
    
    for i = 1:3
        curIds = ceil(borders / prod(sizePadded(1:(i-1))));
        curIds = 1 + mod(curIds - 1, sizePadded(i));
        
        curMask = curIds <= (sizePadded(i) - sphereRadVx(i));
        curMask = curMask & ((1 + sphereRadVx(i)) <= curIds);
        borders = borders(curMask);
    end
    
    borders = Util.indToSubMat(sizePadded, borders);
    borders = Util.sub2ind(sizeRaw, borders - sphereRadVx);
    
    mask = seg(borders);
    mask = mask ~= 1 & mask ~= 2;
    borders = borders(mask);
    
    mask = false(size(seg));
    mask(borders) = true;
    
    isoSurf = isosurface(mask, 0.5);
    isoSurf.faces = reshape(isoSurf.faces, [], 3);
    isoSurf.vertices = reshape(isoSurf.vertices, [], 3);
    
    isoSurf.vertices = isoSurf.vertices(:, [2, 1, 3]);
    isoSurf.vertices = isoSurf.vertices + box(:, 1)';
    
    info = opts.info;
    Util.save(outFile, info, isoSurf);
end
