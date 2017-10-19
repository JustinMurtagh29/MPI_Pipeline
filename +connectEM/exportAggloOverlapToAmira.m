function exportAggloOverlapToAmira( ...
        param, agglos, skelToAgglos, outDir, varargin)
    % exportAggloOverlapToAmira(param, agglos, skelToAgglos, outDir, ...)
    %   Exports all agglomerates which with a given ground truth skeleton
    %   to PLY files for Amira.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % calculate all isosurfaces in parallel
    cellIds = repelem( ...
        (1:numel(skelToAgglos))', ...
        cellfun(@numel, skelToAgglos));
    allAgglos = agglos(cell2mat(skelToAgglos));
    
    allAggloCount = numel(allAgglos);
    isoSurfs = cell(allAggloCount, 1);
    
    for curIdx = 1:allAggloCount
        isoSurfs{curIdx} = Visualization.buildIsoSurface( ...
            param, allAgglos{curIdx}, varargin{:});
    end
    
    % write PLY files
    for curCellIdx = 1:max(cellIds)
        curMask = (cellIds == curCellIdx);
        curIsoSurfs = isoSurfs(curMask);
        
        curOutFile = sprintf('cell-%d.ply', curCellIdx);
        fprintf('Writing %s ... ', curOutFile);
        
        Visualization.exportIsoSurfaceToAmira( ...
            param, curIsoSurfs, fullfile(outDir, curOutFile));
        fprintf('done!\n');
    end
end