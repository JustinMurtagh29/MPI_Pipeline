function dendT = splitWholeCellInputs(wcT, splitNmlT)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
   [~, curIds] = ismember(splitNmlT.cellId, wcT.id);
    dendT = cell(size(splitNmlT.cellId));
   
    for curIdx = 1:height(wcT)
        dendT{curIdx} = forCell( ...
            wcT(curIds(curIdx), :), ...
            splitNmlT(curIdx, :));
    end
    
    dendT = cat(1, dendT{:});
end

function dendT = forCell(wcT, splitNmlT)
    somaNodeIds = find(ismember( ...
        wcT.agglo.nodes(:, 4), ...
        wcT.somaAgglo.nodes(:, 4)));

    dendNodeIds = cellfun( ...
        @(ids) setdiff(ids, somaNodeIds), ...
        splitNmlT.dendNodes{1}, 'UniformOutput', false);
    dendNodeIds(cellfun(@isempty, dendNodeIds)) = [];
    
    dendT = wcT(ones(numel(dendNodeIds), 1), :);
    
    dendT.dendId(:) = 1:numel(dendNodeIds);
    dendT = dendT(:, [1, end, 2:(end - 1)]);
    
    for curId = reshape(dendT.dendId, 1, [])
        curNodeIds = dendNodeIds{curId};
        curNodeIds = setdiff(curNodeIds, somaNodeIds);
        
        dendT.title{curId} = sprintf( ...
            '%s, dendrite %d', dendT.title{curId}, curId);
        dendT.tag{curId} = sprintf( ...
            '%s_dendrite-%d', dendT.tag{curId}, curId);
        
        curAgglo = dendT.agglo(curId);
        curAgglo.nodes = curAgglo.nodes(curNodeIds, :);
       [~, curAgglo.edges] = ismember(curAgglo.edges, curNodeIds);
        curAgglo.edges(~all(curAgglo.edges, 2), :) = [];
        dendT.agglo(curId) = curAgglo;
        
        dendT.nodeDists{curId} = ...
            dendT.nodeDists{curId}(curNodeIds);
        
        curSynapses = dendT.synapses{curId};
       [~, curSynapses.nodeId] = ...
            ismember(curSynapses.nodeId, curNodeIds);
        curSynapses(~curSynapses.nodeId, :) = [];
        dendT.synapses{curId} = curSynapses;
    end
end
