function spineLengths = ...
        calculateSpineLengths(param, trunks, dendrites, spineHeads)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    maxSegId = Seg.Global.getMaxSegId(param);
    trunkLUT = logical(Agglo.buildLUT(maxSegId, trunks));
    
    % Find dendrite to spine head
    attached = Agglo.fromSuperAgglo(dendrites);
    attached = Agglo.buildLUT(maxSegId, attached);
    
    attached = cellfun(@(ids) mode( ...
        nonzeros(attached(ids))), spineHeads);
    attached(isnan(attached)) = 0;
    
   [uniAggloIds, ~, uniSpineIds] = unique(attached);
   
    if ~uniAggloIds(1)
        uniAggloIds(1) = [];
        uniSpineIds = uniSpineIds - 1;
    end
    
    spineLengths = nan(size(attached));
    for curAggloIdx = 1:numel(uniAggloIds)
        curAggloId = uniAggloIds(curAggloIdx);
        curSpineMask = uniSpineIds == curAggloIdx;
        
        curDendrite = dendrites(curAggloId);
        curSpineHeads = spineHeads(curSpineMask);
        
        spineLengths(curSpineMask) = forAgglo( ...
            param, trunkLUT, curDendrite, curSpineHeads);
    end
end

function spineLengths = forAgglo(param, trunkLUT, dendrite, spineHeads)
    spineLengths = nan(size(spineHeads));
    
    % Project spine heads onto nodes
    spineHeadIds = repelem( ...
        transpose(1:numel(spineHeads)), ...
        cellfun(@numel, spineHeads));
   [~, spineHeadNodeIds] = ismember( ...
        cell2mat(spineHeads), dendrite.nodes(:, 4));
    
    spineHeadIds(~spineHeadNodeIds) = [];
    spineHeadNodeIds(~spineHeadNodeIds) = [];
    
    % Find trunk nodes
    trunkNodeIds = ~isnan(dendrite.nodes(:, 4));
    trunkNodeIds(trunkNodeIds) = ...
        trunkLUT(dendrite.nodes(trunkNodeIds, 4));
    trunkNodeIds = find(trunkNodeIds);
    
    % No trunk nodes? Return NaN.
    if isempty(trunkNodeIds); return; end

    edgeLengths = ...
        dendrite.nodes(dendrite.edges(:, 1), 1:3) ...
      - dendrite.nodes(dendrite.edges(:, 2), 1:3);
    edgeLengths = param.raw.voxelSize .* edgeLengths;
    edgeLengths = sqrt(sum(edgeLengths .* edgeLengths, 2));
    
    spineLengths = graph( ...
        dendrite.edges(:, 1), dendrite.edges(:, 2), ...
        edgeLengths, 'OmitSelfLoops');
    spineLengths = distances( ...
        spineLengths, trunkNodeIds, spineHeadNodeIds);
    
    spineLengths = min(spineLengths, [], 1);
    spineLengths = accumarray(spineHeadIds(:), spineLengths(:), [], @min);
end
