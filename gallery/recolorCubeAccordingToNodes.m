function recolorCubeAccordingToNodes(segParam, segOutParam, skeletonFile, cubeCoord)
    % Recolor a single KNOSSOS cube with overlap and cut out right size

    load(skeletonFile);
    % Read oversegmented KNOSSOS cube
    cube = readKnossosCube(segParam.root, segParam.prefix, cubeCoord, 'uint16=>uint16', '', 'raw', 256);
    % Determine one coordinate of KNOSSOS cube (in order to calculate skeleton node positions in local cube) 
    oneOfCube = cubeCoord * 128 + 1 - 64;
    % Get maximal segmentation ID present in new segmentation (number skeletons + 1 (merger))
    maxID = length(nodes);
    segIds = cell(maxID+1,1); % one bigger for merger (see close to end of function)
    segIdsCount = cell(maxID,1);
    for i=1:length(nodes)
        % Calculate nodes in local coordinates of cubes
        rel_coords = bsxfun(@minus, nodes{i}, oneOfCube-1);
        % Exclude all nodes outside of current cube with overlap
        rel_coords(any(rel_coords < 1,2),:) = [];
        rel_coords(any(rel_coords > 256,2),:) = [];
        % Collect segmentation IDs for all objects
        if ~isempty(rel_coords)
            idx = sub2ind(size(cube), rel_coords(:,1), rel_coords(:,2), rel_coords(:,3));
            nodeSegIds = cube(idx);
            nodeSegIds(nodeSegIds == 0) = [];
            segIds{i} = unique(nodeSegIds);
            for j=1:length(segIds{i})
                segIdsCount{i}(j) = sum(nodeSegIds == segIds{i}(j)); 
            end
        end
    end
    % All IDs & then all with more than one skeleton 'hit'
    emptyIdx = cellfun(@isempty, segIds);
    allIds = cat(1,segIds{~emptyIdx});
    allUniqueIds = unique(allIds);
    n = histc(allIds, allUniqueIds);
    mergedIds = allUniqueIds(n > 1);
    % Decide to which skeleton we assign each of those (winner-take-all)
    for i=1:length(mergedIds)
        % Determine which skeletons are involved (and position in segIDs(Count))
        [segIdx nodeCount] = cellfun(@(x,y)findPositionInCellOfVector(x, y, mergedIds(i)), segIds(1:end-1), segIdsCount);
        temp = segIdx == 0;
        segIdx(temp) = [];
        nodeCount(temp) = [];
        skelIdx = 1:length(segIdsCount);
        skelIdx(temp) = [];
        % Remove Id and count from arrays first
        for j=1:length(skelIdx) 
            segIds{skelIdx(j)}(segIdx(j)) = [];
            segIdsCount{skelIdx(j)}(segIdx(j)) = [];
        end
        % If there is a winner, add him back
        [maxVal, maxIdx] = max(nodeCount);
        if sum(nodeCount == maxVal) == 1
            segIds{skelIdx(maxIdx)}(end+1) = mergedIds(i);
            segIdsCount{skelIdx(maxIdx)}(end+1) = maxVal;
        end
    end
    % Finally: recolor cube according to skeletons
    cubeOut = zeros(size(cube), 'uint16');
    for j=1:length(segIds)
        cubeOut(ismember(cube, segIds{j})) = j;
    end
    % Drop overlap
    cubeOut = cubeOut(65:end-64,65:end-64,65:end-64);
    % Write consolidated results to KNOSSOS cube
    writeKnossosCube(segOutParam.root, segOutParam.prefix, cubeCoord, cubeOut, 'uint16');

end

function [segId, nodeCount] = findPositionInCellOfVector(x, y, id) 
% Seems like a quick fix, should probably be improved
    temp = x == id;
    segId = find(temp);
    nodeCount = y(temp);
    segId(isempty(segId)) = 0;
    nodeCount(isempty(nodeCount)) = 0;
end

