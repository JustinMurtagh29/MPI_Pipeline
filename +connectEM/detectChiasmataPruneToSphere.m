function [thisNodes, thisEdges, nodeIds, thisDistSq] = ...
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
    nodeDistSq = pdist2( ...
        nodes(sphereNodeIds, :), nodes(i, :), 'squaredeuclidean');
    thisMask = (nodeDistSq < (p.sphereRadiusOuter ^ 2));

    sphereNodeIds = sphereNodeIds(thisMask);
    nodeDistSq = nodeDistSq(thisMask);
    clear thisMask;

    % Build subgraph
    nodes = nodes(sphereNodeIds, :);
   [~, edges] = ismember(edges, sphereNodeIds);
    edges = edges(all(edges, 2), :);
   [~, i] = ismember(i, sphereNodeIds);
    assert(i > 0);

    %% Cut out `sphereRadiusInner`
    thisOutMask = (nodeDistSq > (p.sphereRadiusInner ^ 2));

   [~, lut] = Graph.findConnectedComponents( ...
        edges(~all(thisOutMask(edges), 2), :), false);
    thisOutMask(lut ~= lut(i)) = true;
    clear lut;

   [~, lut] = Graph.findConnectedComponents(edges, false);
    thisOutMask(lut ~= lut(i)) = false;
    clear lut;

    %%
    edgeMask = all(thisOutMask(edges), 2);

    nodeIds = unique(edges(edgeMask, :));
    nodeIds = reshape(nodeIds, [], 1);

    thisNodes = nodes(nodeIds, :);
    thisDistSq = nodeDistSq(nodeIds);
    thisEdges = edges(edgeMask, :);

    [~, thisEdges] = ismember(thisEdges, nodeIds);
    assert(all(thisEdges(:)));

    nodeIds = sphereNodeIds(nodeIds);
end
