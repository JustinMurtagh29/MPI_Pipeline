function buildForEdgesInBox(param, dataSet, outDir, box)
    % buildForEdgesInBox(param, dataSet, outDir, box)
    %   Builds NML files to query all edges in a given cube.
    %
    % param
    %   Parameters returned by setParameterSettings
    %
    % dataSet
    %   Name of dataSet (e.g. 'ex145'). This is used to load
    %   the parameters going into the NML header.
    %
    % outDir
    %   Path to directory in which the NML files are placed.
    %
    % box
    %   3x2 matrix. Bounding box in which the edges will be
    %   queried.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % config
    batchSize = 100;

    % check input argument
    assert(all(size(box) == [3, 2]));
    
    % get all edges
    edges = Query.getEdgesInBox(param, box);
    
    edgeCount = size(edges, 1);
    batchCount = ceil(edgeCount / batchSize);
    
    % reshuffle order of edges
    % we don't a region to be dominated by one HiWi
    randIds = randperm(edgeCount);
    edges = edges(randIds, :);
    
    % load segment-to-point mapping
    rootDir = param.saveFolder;
    points = load([rootDir, 'segToPointMap.mat']);
    points = points.segToPointMap;
    
    % build all files
    for curIdx = 1:batchCount
        curStartIdx = ...
            1 + batchSize * (curIdx - 1);
        curEndIdx = min(edgeCount, ...
            curStartIdx + batchSize - 1);
        curBatchSize = ...
            curEndIdx - curStartIdx + 1;
        
        % build skeleton
        curEdges = ...
            edges(curStartIdx:curEndIdx, :);
        
        curNodeCoords = ...
            nan(curBatchSize, 2, 3);
        curNodeCoords(:, 1, :) = ...
            points(curEdges(:, 1), :);
        curNodeCoords(:, 2, :) = ...
            points(curEdges(:, 2), :);
        
        curSkel = Query.buildForEdges( ...
            dataSet, curEdges, curNodeCoords);
        
        % build path to output file
        curFileName = ['batch-', num2str(curIdx), '.nml'];
        curFilePath = fullfile(outDir, curFileName);
        
        curSkel.write(curFilePath);
    end
end