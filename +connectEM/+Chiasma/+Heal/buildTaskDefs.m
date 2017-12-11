function [taskDef, endings] = buildTaskDefs(param, endings, axons)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
	taskDef = table;

    numEndings = size(endings, 1);
    taskDef.position = nan(numEndings, 3);
    taskDef.direction = nan(numEndings, 3);
    taskDef.rotation = nan(numEndings, 3);
    
    for curIdx = 1:numEndings
       [taskDef.position(curIdx, :), ...
        taskDef.direction(curIdx, :), ...
        taskDef.rotation(curIdx, :)] = ...
            forEnding( ...
                param.raw.voxelSize, ...
                axons(endings.aggloId(curIdx)), ...
                endings.nodeId(curIdx));
    end
    
    % remove cases without direction
    endings(~any(taskDef.direction, 2), :) = [];
    taskDef(~any(taskDef.direction, 2), :) = [];
end

function [pos, dir, rot] = forEnding(voxelSize, axon, nodeId)
    %% position
    pos = axon.nodes(nodeId, 1:3);
    
    %% direction
    neighNodeIds = getNeighbours(axon, nodeId);
    
    if numel(neighNodeIds) > 1
        srcNodeIds = neighNodeIds;
    else
        assert(isscalar(neighNodeIds));
        neighNeighNodeIds = getNeighbours(axon, neighNodeIds);
        neighNeighNodeIds = setdiff(neighNeighNodeIds, nodeId);
        
        if isscalar(neighNeighNodeIds)
            srcNodeIds = neighNeighNodeIds;
        else
            srcNodeIds = neighNodeIds;
        end
    end
    
    dir = ...
        axon.nodes(nodeId, 1:3) - ...
        axon.nodes(srcNodeIds, 1:3);
    dir = sum(dir, 1);
    
    if ~any(dir)
        rot = nan(1, 3);
        return;
    end
    
    %% rotation
    rot = dir .* voxelSize;
    rot = rot ./ norm(rot);
    
    % to Euler angles
    rot = connectEM.diffToEulerAngle(rot);
    rot = round(reshape(rot, 1, []));
end

function neighs = getNeighbours(axon, nodeId)
    mask = any(axon.edges == nodeId, 2);
    neighs = setdiff(axon.edges(mask, :), nodeId);
end