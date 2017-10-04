function [isIntersection,nrExits] = detectChiasmataSub(startidx, outputFolder)
    load([outputFolder  'prep.mat']);
    isIntersection = false(size(nodes,1),1);
    nrExits = zeros(size(nodes,1),1);
    for i=startidx:5000:size(nodes,1)
        curNrExits = connectEM.detectChiasmataNodes( ...
            nodes, edges, ones(size(edges,1), 1), p, i);
        
        if curNrExits > 3
            isIntersection(i) = true;
            nrExits(i) = curNrExits;
        end
    end
    save([outputFolder 'temp_' num2str(startidx)],'nrExits', 'isIntersection','-v7.3');
end
