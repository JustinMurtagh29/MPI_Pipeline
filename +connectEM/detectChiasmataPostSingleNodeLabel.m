function [output, queryIdx] = detectChiasmataPostSingleNodeLabel( ...
        p, nodes, nodesV, edges, isIntersection, nrExits)
    
   [cc, centerOfCC] = ...
        connectEM.detectChiasmataNodesCluster(p, nodes, isIntersection);
    
    % Find out where to query for each CC
    queryIdx = cell(length(cc),1);
    pos = cell(length(cc),1);
    dir = cell(length(cc),1);
    
    for i = 1:length(cc)
        [~, pos{i}, dir{i}, queryIdx{i}] = ...
            connectEM.detectChiasmataNodes( ...
                p, nodes, edges, cc{i}(centerOfCC(i)));
    end

    % Create an output structure
    output.nodes = nodesV;
    output.edges = edges;
    output.prob = [];
    output.isIntersection = isIntersection;
    output.nrExits = nrExits;
    output.ccNodeIdx = cc;
    output.ccCenterIdx = cellfun(@(x,y)x(y), cc, num2cell(centerOfCC));
    output.queryIdx = queryIdx;
    output.position = pos;
    output.direction = dir;
end
