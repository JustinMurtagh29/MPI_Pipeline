function exportAggloOverlapToAmira(param, aggloOverlaps, outDir, varargin)
    % exportAggloOverlapToAmira(param, aggloOverlaps, outDir)
    %   Exports all agglomerates which with a given ground truth skeleton
    %   to PLY files for Amira.
    %
    % aggloOverlaps
    %   Nested cell arrays where aggloOverlaps{cellIdx}{aggloIdx} contains
    %   the segment IDs in agglomerate `aggloIdx` which overlaps with cell
    %   `cellIdx`.
    %
    % outDir
    %   Output directory. One PLY file per cell in `aggloOverlaps` will be
    %   written. This means that one PLY file might contain multiple
    %   meshes.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    assert(iscell(aggloOverlaps));
    assert(all(@iscell, aggloOverlaps));
    
    % calculate all isosurfaces in parallel
    cellIds = repelem( ...
        (1:numel(aggloOverlaps))', ...
        cellfun(@numel, aggloOverlaps));
    agglos = cat(1, aggloOverlaps{:});
    
    isoSurfs = buildIsoSurfaces( ...
        param, agglos, outDir, varargin{:});
    
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

function isoSurfs = buildIsoSurfaces(param, agglos, varargin)
    cluster = Cluster.getCluster( ...
        '-l h_vmem=12G', ...
        '-l s_rt=11:59:00', ...
        '-l h_rt=12:00:00');
    job = Cluster.startJob( ...
        @Visualization.buildIsoSurface, ...
        cellfun(@(segIds) {{segIds}}, agglos), ...
        'sharedInputs', [{param}, varargin], ...
        'sharedInputLocations', [1, 3:(3 + numel(varargin))], ...
        'numOutputs', 1, 'cluster', cluster, 'name', mfilename());
    
    Cluster.waitForJob(job);
    isoSurfs = fetchOutputs(job);
end