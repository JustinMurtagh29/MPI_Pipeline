function [isIntersection,nrExits] = detectChiasmataSub(startidx, outputFolder)
    load([outputFolder  'prep.mat']);
    isIntersection = false(size(nodes,1),1);
    nrExits = zeros(size(nodes,1),1);
    for i=startidx:5000:size(nodes,1)
        [thisNodes, thisEdges] = connectEM.detectChiasmataPruneToSphere( ...
            nodes, edges, ones(size(edges,1), 1), p, i);
        
        C = Graph.findConnectedComponents(thisEdges, false);
        curNrExits = sum(cellfun(@(idx) ...
            max(pdist2(thisNodes(idx, :), nodes(i, :))) > 3000, C));
        
        if curNrExits > 3
            isIntersection(i) = true;
            nrExits(i) = curNrExits;
        end
    end
    save([outputFolder 'temp_' num2str(startidx)],'nrExits', 'isIntersection','-v7.3');
end
