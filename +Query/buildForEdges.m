function skel = buildForEdges(dataSet, edges, nodeCoords)
    % skel = buildForEdges(edges, nodeCoords)
    %   Build a skeleton from the given list of edges and
    %   node coordinates.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    edgeCount = size(edges, 1);
    
    % build colors
    colorVals = ...
        (1:edgeCount) ...
      * (1 - 2 / (1 + sqrt(5)));
    colorVals = mod(colorVals, 1);
    
    colors = [ ...
        colorVals(:), ...
        ones(edgeCount, 1), ...
        ones(edgeCount, 1)];
    
    % build RGBA values
    colors = hsv2rgb(colors);
    colors = [colors, ones(edgeCount, 1)];
    
    % sanity check
    assert(edgeCount == size(edges, 1));
    
    % build skeleton
    skel = skeleton();
    
    % build all trees
    for curIdx = 1:edgeCount
        curEdge = edges(curIdx, :);
        curColor = colors(curIdx, :);
        
        curCoords = squeeze( ...
            nodeCoords(curIdx, :, :));
        assert(all(size(curCoords) == [2, 3]));
        
        % build name
        curName = sprintf( ...
            'Edge %d-%d', curEdge(1), curEdge(2));
        
        % add tree to skeleton
        skel = skel.addTree( ...
            curName, curCoords, [1, 2], curColor);
    end
    
    % set active node and zoom
    skel.parameters.activeNode.id = '1';
    skel.parameters.zoomLevel.zoom = '0.5604482797441586';
    
    % set edit position
    skel.parameters.editPosition.x = num2str(nodeCoords(1, 1));
    skel.parameters.editPosition.y = num2str(nodeCoords(1, 2));
    skel.parameters.editPosition.z = num2str(nodeCoords(1, 3));
        
    % set parameter
    skel = Skeleton.setParams4Dataset(skel, dataSet);
end