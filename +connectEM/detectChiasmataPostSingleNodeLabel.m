function [output, queryIdx] = ...
        detectChiasmataPostSingleNodeLabel( ...
            p, nodes, nodesV, edges, isIntersection, nrExits)
    % Written by
    %   Manuel Berning <manuel.berning@brain.mpg.de>
    %   Kevin Boergens <kevin.boergens@brain.mpg.de>
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % Cluster nodes
   [clusters, centerIds] = ...
        connectEM.Chiasma.Detect.clusterNodes( ...
            nodes, isIntersection, p.clusterSize);
    
    % Generate queries
   [~, pos, dir, queryIdx] = arrayfun(@(i) ...
       connectEM.detectChiasmataNodes(p, nodes, edges, id), ...
       centerIds, 'UniformOutput', false);

    % Create an output structure
    output.nodes = nodesV;
    output.edges = edges;
    output.prob = zeros(0, 1);
    output.isIntersection = isIntersection;
    output.nrExits = nrExits;
    output.ccNodeIdx = clusters;
    output.ccCenterIdx = centerIds;
    output.queryIdx = queryIdx;
    output.position = pos;
    output.direction = dir;
end
