function asiAreas = buildAxonSpineInterfaceAreas( ...
        param, graph, axons, shAgglos, synapses, synT)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    maxSegId = Seg.Global.getMaxSegId(param);
    shLUT = Agglo.buildLUT(maxSegId, shAgglos);
    
    synT.shId = cellfun( ...
        @(segIds) max(shLUT(segIds)), ...
        synapses.postsynId(synT.id));
    uniShIds = setdiff(synT.shId, 0);
   [~, synT.shId] = ismember(synT.shId, uniShIds);
	
    %% Build list of outward edges for each spine head.
    shAgglos = shAgglos(uniShIds);
    shLUT = Agglo.buildLUT(maxSegId, shAgglos);
    
    shEdgesIds = [ ...
        shLUT(graph.edges(:, 1)), reshape(1:size(graph., 1), [], 1); ...
        shLUT(graph.edges(:, 2)), reshape(1:size(graph, 1), [], 1)];
    shEdgesIds(~shEdgesIds(:, 1), :) = [];
    
    shEdgesIds = accumarray( ...
        shEdgesIds(:, 1), shEdgesIds(:, 2), ...
        [], @(edgeIds) {sort(edgeIds)});
    
    %% Calculate axon spine interface areas
    asiAreas = nan(size(synT.id));
    for curIdx = 1:size(synT, 1)
        if ~synT.isSpine(curIdx)
            continue;
        end
        
        curShEdgeIds = shEdgesIds{synT.shId(curIdx)};
        curPreSegIds = axons{synT.preAggloId(curIdx)};
        
        curAsiArea = graph.edges(curShEdgeIds, :);
        curAsiArea = any(ismember(curAsiArea, curPreSegIds), 2);
        curAsiArea = sum(graph.borderArea(curShEdgeIds(curAsiArea)));
        
        asiAreas(curIdx) = curAsiArea;
    end
end