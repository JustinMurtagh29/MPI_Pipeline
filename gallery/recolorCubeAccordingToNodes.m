function recolorCubeAccordingToNodes(segParam, segOutParam, skeletonFile, cubeCoord)
    % Recolor a single KNOSSOS cube with overlap and cut out right size

    load(skeletonFile);
    % Read oversegmented KNOSSOS cube
    cube = readKnossosCube(segParam.root, segParam.prefix, cubeCoord, 'uint16=>uint16', '', 'raw', 256);
    % Determine one coordinate of KNOSSOS cube (with overlap for reading and without for writing) 
    oneOfCube = cubeCoord * 128 + 1;
    % Get maximal segmentation ID present in new segmentation (number skeletons + 1 (merger))
    maxID = length(nodes) + 1;
    segIds = cell(maxID,1);
    for i=1:length(nodes)
        % Calculate nodes in local coordinates of cubes
        rel_coords = bsxfun(@minus, nodes{i} + 1, oneOfCube);
        % Exclude all nodes outside of current cube with overlap
        rel_coords(any(rel_coords < 1,2),:) = [];
        rel_coords(any(rel_coords > 256,2),:) = [];
        % Collect segmentation IDs for all objects
        if ~isempty(rel_coords)
            idx = sub2ind(size(cube), rel_coords(:,1), rel_coords(:,2), rel_coords(:,3));
            segIds{i} = unique(cube(idx));
            segIds{i}(segIds{i} == 0) = [];
        end
    end
    % Find merger
    for i=1:length(segIds)
        for j=(i+1):length(segIds)
            merger{i,j} = intersect(segIds{i},segIds{j});
        end
    end
    mergerFound = ~cellfun(@isempty, merger);
    % Set all mergers found to next free id (number skeletons + 1)
    [row, col] = find(mergerFound == 1);
    for i=1:length(row)
        ids = merger{row(i),col(i)};
        segIds{row(i)}(ismember(segIds{row(i)}, ids)) = [];
        segIds{col(i)}(ismember(segIds{col(i)}, ids)) = [];
        segIds{end} = [segIds{end}; ids];
    end
    % Recolor cube according to skeletons
    cubeOut = zeros(size(cube), 'uint16');
    for j=1:length(segIds)
        cubeOut(ismember(cube, segIds{j})) = j;
    end 
    % Drop overlap
    cube = cube(65:end-64,65:end-64,65:end-64);
    % Write consolidated results to KNOSSOS cube
    writeKnossosCube(segOutParam.root, segOutParam.prefix, cubeCoord, cubeOut, 'uint16');

end

