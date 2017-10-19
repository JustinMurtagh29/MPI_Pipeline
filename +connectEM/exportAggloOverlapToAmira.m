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
    
    isoSurfs = buildIsoSurfaces( ...
        param, allAgglos, outDir, varargin{:});
    
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

function isoSurfs = buildIsoSurfaces(param, agglos, outDir, varargin)
    tempDir = fullfile(outDir, 'agglomerates');
    mkdir(tempDir);
    
    outFiles = arrayfun(@(idx) ...
        fullfile(tempDir, sprintf('agglo-%d.mat', idx)), ...
        1:numel(agglos), 'UniformOutput', false);
    inputs = cellfun(@(segIds, outFile) ...
        {{segIds, outFile}}, agglos(:), outFiles(:));
    sharedInputs = cat(2, {param}, varargin);
        
    cluster = Cluster.getCluster( ...
        '-l h_vmem=12G', ...
        '-l s_rt=3:59:00', ...
        '-l h_rt=4:00:00');
    
    job = Cluster.startJob( ...
        @jobFunction, inputs, ...
        'sharedInputs', sharedInputs, ...
        'sharedInputsLocation', [1, 3:(3 + numel(varargin))], ...
        'cluster', cluster, 'name', mfilename());
    wait(job);
    
    isoSurfs = cellfun(@load, outFiles(:), 'UniformOutput', false);
end

function jobFunction(varargin)
    isoSurf = Visualization.buildIsoSurface(varargin{1:(end - 1)});
    Util.saveStruct(varargin{end}, isoSurf);
end