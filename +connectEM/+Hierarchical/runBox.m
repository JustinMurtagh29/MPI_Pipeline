function [mergedEdges, mergeScores] = runBox(param, box, varargin)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opt = struct;
    opt.margin = [0; 0; 0];
    opt = Util.modifyStruct(opt, varargin{:});
    
    boxLarge = box + [-1, +1] .* opt.margin(:);
    boxLarge = max(boxLarge, param.bbox(:, 1));
    boxLarge = min(boxLarge, param.bbox(:, 2));
    
    seg = loadSegDataGlobal(param.seg, boxLarge);
    coreSegIds = findCoreSegIds(box, boxLarge, seg);
    
    maxSegId = Seg.Global.getMaxSegId(param);
    corrEdges = loadCorrespondences(param, boxLarge);
    corrLUT = buildCorrespondenceLUT(maxSegId, corrEdges);
    seg(seg ~= 0) = corrLUT(seg(seg ~= 0));
    
   [edges, borders] = SynEM.Svg.findEdgesAndBorders(seg);
    borders = reshape({borders.PixelIdxList}, [], 1);
    class = loadClassData(param.class, boxLarge);
    
   [mergedBorders, mergeScores] = ...
        connectEM.Hierarchical.run(class, edges, borders);
    
    mask = any(ismember(edges, coreSegIds), 2);
    mask = cellfun(@(ids) any(mask(ids)), mergedBorders);

    mergedEdges = [{corrEdges}; cellfun( ...
        @(ids) unique(edges(ids, :), 'rows'), ...
        mergedBorders(mask), 'UniformOutput', false)];
    mergeScores = [inf; mergeScores(mask)];
end

function coreSegIds = findCoreSegIds(box, boxLarge, seg)
    if ~isequal(box, boxLarge)
        off = box - boxLarge;
        seg = seg( ...
            (1 + off(1, 1)):(end + off(1, 2)), ...
            (1 + off(2, 1)):(end + off(2, 2)), ...
            (1 + off(3, 1)):(end + off(3, 2)));
    end
    
    coreSegIds = setdiff(seg, 0);
    coreSegIds = reshape(coreSegIds, [], 1);
end

function lut = buildCorrespondenceLUT(maxSegId, corr)
    lut = zeros(maxSegId, 1, 'like', corr);
    lut(1:end) = 1:numel(lut);
    
    % NOTE(amotta): This assumes that `corr` are sorted both column- and
    % row-wise. This should be true for the result of `loadCorrespondences`.
    for i = 1:size(corr, 1); lut(corr(i, 2)) = lut(corr(i, 1)); end
end

function corrEdges = loadCorrespondences(param, box)
    corrDir = param.correspondence.saveFolder;
    
    tileIds = box - param.bbox(:, 1) + 1;
    tileIds = ceil(tileIds ./ param.tileSize);
    
    tileCount = tileIds(:, 2) - tileIds(:, 1) + 1;
    corrEdges = repmat({zeros(0, 2, 'uint32')}, [tileCount', 3]);
    
    for curDim = 1:3
        for curX = 1:(tileCount(1) - (curDim == 1))
            for curY = 1:(tileCount(2) - (curDim == 2))
                for curZ = 1:(tileCount(3) - (curDim == 3))
                   [curIdsOne, curIdsTwo] = deal([ ...
                        tileIds(1) + curX - 1, ...
                        tileIds(2) + curY - 1, ...
                        tileIds(3) + curZ - 1]);
                    curIdsTwo(curDim) = curIdsTwo(curDim) + 1;
                    
                    curCorr = sprintf( ...
                        '%s%sglobal.mat', ...
                        num2str(curIdsOne, '%02d'), ...
                        num2str(curIdsTwo, '%02d'));
                    curCorr = fullfile(corrDir, curCorr);
                    curCorr = load(curCorr);
                    
                    % NOTE(amotta): Only use correspondences where the
                    % contact area is at least 10 vx and the contact makes
                    % up at least 90 % of the smaller segment-border areas.
                    %
                    % See +connectEM/collectGlobalCorrespondences.m in
                    % commit 20ed1a080f50397f32c19f77d57eb9c8cba92c57.
                   [~, curMask] = ismember( ...
                       curCorr.uniqueCorrespondences, ...
                       curCorr.uniqueSegments);
                    curMask = curCorr.countsS(curMask);
                    curMask = reshape(curMask, [], 2);
                    curMask = min(curMask, [], 2);

                    curMask = ...
                        curCorr.countsC > 10 ...
                      & curCorr.countsC ./ curMask > 0.9;
                  
                    curCorr = curCorr.uniqueCorrespondences(curMask, :);
                    corrEdges{curX, curY, curZ, curDim} = curCorr;
                end
            end
        end
    end
    
    corrEdges = cell2mat(corrEdges(:));
    corrEdges = unique(sort(corrEdges, 2), 'rows');
end

