function edges = getEdgesInBox(param, box)
    % edges = getEdgesInBox(param, box)
    %   Returns a list of all edges that are fully contained
    %   in the user-specified bounding box.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    rootDir = param.saveFolder;
    
    % load segmentation
    seg = loadSegDataGlobal( ...
        param.seg.root, param.seg.prefix, box);
    
    % get segment ids
    segIds = unique(seg(:));
    segIds = segIds(segIds > 0);
    
    % load all edges
    graph = load([rootDir, 'graph.mat'], 'edges');
    
    % find edges in box
    keepEdgesMask = all( ...
        ismember(graph.edges, segIds), 2);
    edges = graph.edges(keepEdgesMask, :);
    edges = unique(edges, 'rows');
end