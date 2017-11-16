function [thisNodes, thisEdges, nodeIds, thisDist] = ...
        detectChiasmataPruneToSphere(p, nodes, edges, i)
    % Written by
    %   Kevin Boergens <kevin.boergens@brain.mpg.de>
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    assert(isscalar(i));
    
    %% Restrict to `sphereRadiusOuter`
    % First pass, using bounding box
    sphereNodeIds = find( ...
        all(bsxfun(@ge, nodes, nodes(i, :) - p.sphereRadiusOuter), 2) ...
      & all(bsxfun(@le, nodes, nodes(i, :) + p.sphereRadiusOuter), 2));

    % Second pass, using `pdist2`
    nodeDist = pdist2(nodes(sphereNodeIds, :), nodes(i, :));
    thisMask = (nodeDist < p.sphereRadiusOuter);

    sphereNodeIds = sphereNodeIds(thisMask);
    nodeDist = nodeDist(thisMask);
    clear thisMask;

    % Build subgraph
    nodes = nodes(sphereNodeIds, :);
   [~, edges] = ismember(edges, sphereNodeIds);
    edges = edges(all(edges, 2), :);
   [~, i] = ismember(i, sphereNodeIds);
    assert(i > 0);

    %% Cut out `sphereRadiusInner`
    thisOutMask = (nodeDist > p.sphereRadiusInner);
    thisOutMask = reshape(thisOutMask, 1, []);

   [~, lut] = Graph.findConnectedComponents( ...
        edges(~all(thisOutMask(edges), 2), :), ...
        false, false, size(nodes, 1));
    thisOutMask(lut ~= lut(i)) = true;
    clear lut;

   [~, lut] = Graph.findConnectedComponents( ...
        edges, false, false, size(nodes, 1));
    thisOutMask(lut ~= lut(i)) = false;
    clear lut;

    %%
    edgeMask = all(thisOutMask(edges), 2);

    nodeIds = unique(edges(edgeMask, :));
    nodeIds = reshape(nodeIds, [], 1);

    thisNodes = nodes(nodeIds, :);
    thisDist = nodeDist(nodeIds);
    thisEdges = edges(edgeMask, :);

    [~, thisEdges] = ismember(thisEdges, nodeIds);
    assert(all(thisEdges(:)));

    nodeIds = sphereNodeIds(nodeIds);
end
