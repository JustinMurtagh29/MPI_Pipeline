function skel = buildSkeleton(agglo, chiasmata, chiasmaIds)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    chiasmaCount = numel(chiasmata.ccCenterIdx);
    nodeCount = size(agglo.nodes, 1);
    
    if ~exist('chiasmaIds', 'var') || isempty(chiasmaIds)
        % by default, show all chiasmata
        chiasmaIds = 1:chiasmaCount;
    end
    
    if ~isfield(agglo, 'solvedChiasma')
        agglo.solvedChiasma = false(nodeCount, 1);
    end
    
    % build comments
    comments = repelem({''}, nodeCount, 1);
    solvedStr = {'unsolved', 'solved'};
    
    for curIdx = reshape(chiasmaIds, 1, [])
        curNodeId = chiasmata.ccCenterIdx(curIdx);
        curSolved = agglo.solvedChiasma(curNodeId);
        
        comments{curNodeId} = sprintf( ...
            'Chiasma %d (%s)', curIdx, solvedStr{1 + curSolved});
        
        comments(chiasmata.queryIdx{curIdx}) = arrayfun( ...
            @(exitIdx) sprintf('Chiasma %d, Exit %d', curIdx, exitIdx), ...
            1:numel(chiasmata.queryIdx{curIdx}), 'UniformOutput', false);
    end
    
    skel = skeleton();
    skel = skel.addTree( ...
        sprintf('Agglomerate (%d chiasmata)', chiasmaCount), ...
        agglo.nodes(:, 1:3), agglo.edges, [], [], comments);
    skel = skel.addBranchpoint( ...
        find(not(cellfun(@isempty, comments)))); %#ok
end