function [spineLengths, dendIds, nodeIds] = ...
        calculateSpineLengths(param, trunks, dendrites, spineHeads)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    spineHeads = reshape(spineHeads, [], 1);
    
    maxSegId = Seg.Global.getMaxSegId(param);
    trunkLUT = logical(Agglo.buildLUT(maxSegId, trunks));
    
    % Find dendrite to spine head
    dendIds = Agglo.fromSuperAgglo(dendrites);
    dendIds = Agglo.buildLUT(maxSegId, dendIds);
    dendIds = cellfun(@(ids) mode(nonzeros(dendIds(ids))), spineHeads);
    dendIds(isnan(dendIds)) = 0;
    
   [uniAggloIds, ~, uniSpineIds] = unique(dendIds);
   
    if ~uniAggloIds(1)
        uniAggloIds(1) = [];
        uniSpineIds = uniSpineIds - 1;
    end
    
    spineLengths = nan(size(spineHeads));
    nodeIds = nan(size(spineHeads));
    for curAggloIdx = 1:numel(uniAggloIds)
        curAggloId = uniAggloIds(curAggloIdx);
        curSpineMask = uniSpineIds == curAggloIdx;
        
        curDendrite = dendrites(curAggloId);
        curSpineHeads = spineHeads(curSpineMask);
        
       [curLengths, curNodeIds] = forAgglo( ...
            param, trunkLUT, curDendrite, curSpineHeads);
        
        spineLengths(curSpineMask) = curLengths;
        nodeIds(curSpineMask) = curNodeIds;
    end
end

function [spineLengths, nodeIds] = forAgglo( ...
        param, trunkLUT, dendrite, spineHeads)
    spineLengths = nan(size(spineHeads));
    nodeIds = nan(size(spineHeads));
    
    flatten = @(v) reshape(v, [], 1);
    
    shT = table;
    shT.id = flatten(repelem( ...
        1:numel(spineHeads), ...
        cellfun(@numel, spineHeads)));
    
    shT.nodeId = cell2mat(spineHeads(:));
    dendriteSegIds = dendrite.nodes(:, 4);
   [~, shT.nodeId] = ismember(shT.nodeId, dendriteSegIds);
    shT = shT(shT.nodeId > 0, :);
    
    % Find trunk nodes
    trunkNodeIds = not( ...
        isnan(dendriteSegIds));
    trunkNodeIds(trunkNodeIds) = ...
        trunkLUT(dendriteSegIds(trunkNodeIds));
    trunkNodeIds = find(trunkNodeIds);
    
    % No trunk nodes? Return NaN.
    if isempty(trunkNodeIds); return; end

    trunkDists = ...
        dendrite.nodes(dendrite.edges(:, 1), 1:3) ...
      - dendrite.nodes(dendrite.edges(:, 2), 1:3);
    trunkDists = param.raw.voxelSize .* trunkDists;
    trunkDists = sqrt(sum(trunkDists .* trunkDists, 2));
    
    trunkDists = graph( ...
        dendrite.edges(:, 1), dendrite.edges(:, 2), ...
        trunkDists, size(dendrite.nodes, 1), 'OmitSelfLoops');
    trunkDists = distances(trunkDists, trunkNodeIds, shT.nodeId);
    
   [shT.length(:), shT.trunkNodeId(:)] = min(trunkDists, [], 1);
    shT.trunkNodeId = trunkNodeIds(shT.trunkNodeId);
    
    shT = sortrows(shT, 'length', 'ascend');
   [shIds, shRows] = unique(shT.id, 'stable');
    spineLengths(shIds) = shT.length(shRows);
    nodeIds(shIds) = shT.trunkNodeId(shRows);
end
