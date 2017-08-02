function pathLength = getPathLengthFromAgglomeration(idsCell, com)
    % Calculate path length of minimal spanning tree of CoM of segments of all cells in idsCell
    comCell = cellfun(@(x)bsxfun(@times, com(x,:), [11.24 11.24 28]), idsCell, 'UniformOutput', false);        
    edgesCell = cellfun(@(x) Graph.getMST(x,2000), comCell, 'UniformOutput', false);  % 2000 nm edge threshold
    pathLength = cellfun(@(x,y)getPathLength(x,y), comCell, edgesCell);
end

function pathLength = getPathLength(com, edges)
    if isempty(edges)
        pathLength = 0;
    else
        edgePos1 = com(edges(:,1),:);
        edgePos2 = com(edges(:,2),:);
        lengths = sqrt(sum((edgePos1-edgePos2).^2,2));
        pathLength = sum(lengths)./1e3;
    end
end


