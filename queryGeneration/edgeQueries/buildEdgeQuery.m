function skel = buildEdgeQuery(com, edge,experimentName)
    % skel = buildEdgeQuery(com, edge)
    %   Build skeleton with query for specified edge.
    %
    % com
    %   N x 3 double matrix. Row i contains the x, y and z
    %   coordinates for the global segment ID i.
    %
    % edge
    %   1 x 2 matrix. Each entry is a global segment ID.
    % 
    % experimentName
    %   String specifying the name of dataset (see switch-case 
    %   in Skeleton.setParams4Dataset in auxiliary methods)
    %   For Example: experimentName = 'ex145'
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    % Modifed by
    %   Sahil Loomba
    
    edgeName = sprintf( ...
        'Edge %d-%d', edge(1), edge(2));
    
    queryNodes = [ ...
        com(edge(1), :); ...
        com(edge(2), :)];
    
    % build skeleton
    skel = skeleton();
    skel = skel.addTree(edgeName, queryNodes);
    
    % set active node and zoom
    skel.parameters.activeNode.id = '1';
    skel.parameters.zoomLevel.zoom = '0.5604482797441586';
    
    % set edit position
    skel.parameters.editPosition.x = num2str(queryNodes(1, 1));
    skel.parameters.editPosition.y = num2str(queryNodes(1, 2));
    skel.parameters.editPosition.z = num2str(queryNodes(1, 3));
    
    % set colors
    skel.colors{1} = [1, 1, 0, 1];
    
    % set parameter
    skel = Skeleton.setParams4Dataset(skel, experimentName);
end
