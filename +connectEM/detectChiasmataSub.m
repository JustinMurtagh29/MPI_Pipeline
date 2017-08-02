function [isIntersection,nrExits] = detectChiasmataSub(startidx)
    load(['/tmpscratch/kboerg/visX17_0/visX17_1/prep.mat']);
    isIntersection = false(size(nodes,1),1);
    nrExits = zeros(size(nodes,1),1);
    for i=startidx:5000:size(nodes,1)
        [thisNodes, thisEdges, thisProb] = connectEM.detectChiasmataPruneToSphere(nodes, edges, ones(size(edges,1),1), p, i);
        C = Graph.findConnectedComponents(thisEdges);
        if length(C) > 3 && sum(cellfun(@(idx)max(pdist2(thisNodes(idx, :), nodes(i,:))) > 4000, C))>3
            isIntersection = true;
            nrExits = length(C);
        else
            isIntersection = false;
            nrExits = 0;
        end
    end
    save([outputFolder 'temp_' num2str(startidx)],'nrExits', 'isIntersection');
end
