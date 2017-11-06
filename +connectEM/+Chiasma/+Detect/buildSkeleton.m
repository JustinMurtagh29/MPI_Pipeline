function skel = buildSkeleton(agglo, chiasmata, chiasmaIds)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    chiasmaCount = numel(chiasmata.ccCenterIdx);
    
    if ~exist('chiasmaIds', 'var') || isempty(chiasmaIds)
        % by default, show all chiasmata
        chiasmaIds = 1:chiasmaCount;
    end
    
    % build comments
    nodeCount = size(agglo.nodes, 1);
    comments = repelem({''}, nodeCount, 1);
    
    for chiIdx = reshape(chiasmaIds, 1, [])
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