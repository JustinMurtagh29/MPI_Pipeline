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
    ballRadiusInNm = 160;
    
    % get parameters
    cubeParam = param.local(cubeIdx);
    cubeDir = cubeParam.saveFolder;
    
    % load edges and borders
    load([cubeDir, 'edges.mat'], 'edges');
    load([cubeDir, 'borders.mat'], 'borders');
    
    % load segmentation data
    seg = loadSegData(param, cubeIdx);
    
    % build ball offsets
    offsetIds = Border.extensionBall( ...
        param, ballRadiusInNm);
    
    % prepare output
    edgeCount = size(edges, 1);
    bordersExt = cell(edgeCount, 3);
    
    doFunc = @(edge, border) ...
        extendBorder(seg, offsetIds, edge, border);
    
    tic;
    for curIdx = 1:edgeCount
        curBorder = borders(curIdx);
        curBorderIds = curBorder.PixelIdxList;
        curBorderIds = curBorderIds(:);
        
        [curIdsOne, curIdsTwo] = ...
            doFunc(edges(curIdx, :), curBorderIds);
        
        % save result
        bordersExt{curIdx, 1} = curBorderIds;
        bordersExt{curIdx, 2} = curIdsOne;
        bordersExt{curIdx, 3} = curIdsTwo;
        
        % show progress
        Util.progressBar(curIdx, edgeCount);
    end
    
    % Prepare result
    outStruct = struct;
    outStruct.borders = borderExt;
    
    outFile = [cubeDir, 'edgesExt.mat'];
    save(outFile, '-struct', outStruct);
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

function data = loadSegData(param, cubeIdx)
    cubeParam = param.local(cubeIdx);
    cubeBox = cubeParam.bboxSmall;
    
    % load data
    data = loadSegDataGlobal( ...
        param.seg.root, param.seg.prefix, cubeBox);
end