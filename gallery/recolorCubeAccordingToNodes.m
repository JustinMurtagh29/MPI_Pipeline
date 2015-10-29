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
    % Find merger
    for i=1:length(segIds)
        for j=(i+1):length(segIds)
            merger{i,j} = intersect(segIds{i},segIds{j});
        end
    end
    mergerFound = ~cellfun(@isempty, merger);
    % Look at each merger and do winner-take-all
    [row, col] = find(mergerFound);
    for i=1:length(row)
        ids = merger{row(i),col(i)};
        for j=1:length(ids)
            idxRow = ismember(segIds{row(i)}, ids(j)); 
            idxCol = ismember(segIds{col(i)}, ids(j));
            countRow = segIdsCount{row(i)}(idxRow);
            countCol = segIdsCount{col(i)}(idxCol);
            if countRow > countCol
                segIds{col(i)}(idxCol) = [];
            elseif countCol > countRow
                segIds{row(i)}(idxRow) = [];
            else
                % If no winner is found, set to maxID + 1 (unresolved mergers)
                segIds{row(i)}(idxRow) = [];
                segIds{col(i)}(idxCol) = [];
                segIds{maxID+1} = [segIds{maxID+1}; ids];
            end
        end
    end
    % Recolor cube according to skeletons
    cubeOut = zeros(size(cube), 'uint16');
    for j=1:length(segIds)
        cubeOut(ismember(cube, segIds{j})) = j;
    end
    % Drop overlap
    cubeOut = cubeOut(65:end-64,65:end-64,65:end-64);
    % Write consolidated results to KNOSSOS cube
    writeKnossosCube(segOutParam.root, segOutParam.prefix, cubeCoord, cubeOut, 'uint16');

end

