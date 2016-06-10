function extend(param, cubeIdx)
    % extend(param, cubeIdx)
    %   Extends all borders in a segmentation cube. For
    %   each edge / border, the following operations are
    %   performed:
    %
    %   1) Lookup voxels in border
    %   2) Find all voxels at most X nm from border
    %   3) Find all voxels that are then contained in the
    %      segments making up the edge
    %
    %   The result is saved to 'bordersExt.mat' and has the
    %   following form:
    %
    %   borders:
    %       Nx3 cell array. Each entry contains linear voxel
    %       indices relative to the small bounding box. Entry
    %       borders{i, j} corresponds to the edges i and
    %
    %       j = 1 for the border
    %       j = 2 for the segment edge(i, 1)
    %       j = 3 for the segment edge(i, 2)
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % config
    ballRadiusInNm = 100;
    
    % get parameters
    cubeParam = param.local(cubeIdx);
    cubeDir = cubeParam.saveFolder;
    cubeBoxSmall = cubeParam.bboxSmall;
    cubeBoxBig = cubeParam.bboxBig;
    
    % load edges and borders
    load([cubeDir, 'edges.mat'], 'edges');
    load([cubeDir, 'borders.mat'], 'borders');
    
    % load segmentation data on the BIG bounding box
    seg = loadSegData(param, cubeIdx);
    
    % transfer linear voxel indices from
    % the SMALL bounding box to the BIG one
    borders = arrayfun(@(b) fixVoxelIds( ...
        cubeBoxSmall, cubeBoxBig, b.PixelIdxList), ...
        borders, 'UniformOutput', false);
    
    % build ball offsets
    [ballOffIds, ballSize] = ...
        Border.extensionBall(param, ballRadiusInNm);
    ballPadding = (ballSize - 1) / 2;
    
    % prepare output
    edgeCount = size(edges, 1);
    bordersExt = cell(edgeCount, 3);
    
    doFunc = @(edge, border) ...
        extendBorder(seg, ballOffIds, edge, border);
    
    tic;
    for curIdx = 1:edgeCount
        curBorderIds = borders{curIdx};
        
        [curIdsOne, curIdsTwo] = ...
            doFunc(edges(curIdx, :), curBorderIds);
        
        % save result
        bordersExt{curIdx, 1} = curBorderIds;
        bordersExt{curIdx, 2} = curIdsOne;
        bordersExt{curIdx, 3} = curIdsTwo;
        
        % show progress
        Util.progressBar(curIdx, edgeCount);
    end
    
    % determine bounding box in which the feature
    % calculate will take place later on
    featBox = cubeBoxSmall;
    featBox(:, 1) = featBox(:, 1) - ballPadding(:);
    featBox(:, 2) = featBox(:, 2) + ballPadding(:);
    
    % fix linear indices to the feature bounding box
    bordersExt = cellfun(@(b) ...
        fixVoxelIds(cubeBoxBig, featBox, b), ...
        bordersExt, 'UniformOutput', false);
    
    % Prepare result
    outStruct = struct;
    outStruct.box = featBox;
    outStruct.borders = bordersExt;
    
    outFile = [cubeDir, 'bordersExt.mat'];
    save(outFile, '-struct', 'outStruct');
end

function [idsOne, idsTwo] = ...
        extendBorder(seg, offsetIds, edge, borderIds)
    % adjust shape
    borderIds = borderIds(:);
    offsetIds = offsetIds(:)';
    
    % find all voxel indices
    allIds = bsxfun( ...
        @plus, borderIds, offsetIds);
    allIds = unique(allIds(:));
    
    % remove voxels outside bounding box
    voxCount = numel(seg);
    allIds = allIds( ...
        (allIds >= 1) ...
      & (allIds <= voxCount));
    
    % find indices in segments
    idsOne = allIds(seg(allIds) == edge(1));
    idsTwo = allIds(seg(allIds) == edge(2));
end

function voxelIds = fixVoxelIds(oldBox, newBox, voxelIds)
    % allocate temporary data
    voxelCount = numel(voxelIds);
    coordMat = nan(voxelCount, 3);
    
    % compute coordinates
    [coordMat(:, 1), coordMat(:, 2), coordMat(:, 3)] = ...
        ind2sub(1 + oldBox(:, 2)' - oldBox(:, 1)', voxelIds(:));
    
    % shift coordinates to new box
    diffVec = newBox(:, 1)' - oldBox(:, 1)';
    coordMat = bsxfun(@minus, coordMat, diffVec);
    
    % to linear indices
    voxelIds = sub2ind( ...
        1 + newBox(:, 2)' - newBox(:, 1)', ...
        coordMat(:, 1), coordMat(:, 2), coordMat(:, 3));
    voxelIds = int32(voxelIds);
end

function seg = loadSegData(param, cubeIdx)
    cubeParam = param.local(cubeIdx);
    cubeDir = cubeParam.saveFolder;
    
    % load global segmentation data
    % for the BIG bounding box
    seg = load([cubeDir, 'segGlobal.mat'], 'seg');
    seg = seg.seg;
end