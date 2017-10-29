function skel = buildSkeletonFromDetection(agglo, chiasmata)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    chiasmaCount = numel(chiasmata.ccCenterIdx);
    nodeCount = size(agglo.nodes, 1);
    
    % build comments
    comments = repelem({''}, nodeCount, 1);
    
    for chiIdx = 1:chiasmaCount
        comments{chiasmata.ccCenterIdx(chiIdx)} = ...
            sprintf('Chiasma %d', chiIdx);
        
        comments(chiasmata.queryIdx{chiIdx}) = arrayfun( ...
            @(exitIdx) sprintf('Chiasma %d, Exit %d', chiIdx, exitIdx), ...
            1:numel(chiasmata.queryIdx{chiIdx}), 'UniformOutput', false);
    end
    
    skel = skeleton();
    skel = skel.addTree( ...
        sprintf('Agglomerate (%d chiasmata)', chiasmaCount), ...
        agglo.nodes(:, 1:3), agglo.edges, [], [], comments);
    skel = skel.addBranchpoint(find(not(cellfun(@isempty, comments)))); %#ok
end