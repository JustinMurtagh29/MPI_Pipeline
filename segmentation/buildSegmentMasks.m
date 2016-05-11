function buildSegmentMasks(cubeParams)
    cubeFolder = cubeParams.saveFolder;
    
    disp('Loading global segment ids');
    segIds = getSegIds(cubeParams);
    sortedIds = sort(segIds);
    
    % find extremal segment ids
    minId = sortedIds(1);
    maxId = sortedIds(end);
    numSegs = maxId - minId + 1;
    
    disp('Loading global segment matrix');
    segGlobalFile = [cubeFolder, 'segGlobal.mat'];
    load(segGlobalFile, 'seg');
    
    % ATTENTION!
    % seg = seg - minId + 1
    % is NOT THE SAME AS
    % seg = seg + 1 - minId = seg - (minId - 1)
    % due to the fact that its uint32
    %
    % make sure the first segment has id one
    seg = seg - (minId - 1);
    
    % calculate region props
    disp('Calculating bounding boxes');
    segProps = regionprops(seg, 'BoundingBox');
    numBoxes = length(segProps);
    
    % sanity check
    if numBoxes ~= numSegs
        error('Sanity check failed!');
    end
    
    % prepare output
    out = struct;
    out.segIds = uint32(sortedIds(:));
    out.boxLocal = uint32(numSegs, 3, 2);
    out.boxGlobal = uint32(numSegs, 3, 2);
    out.masks = cell(numSegs, 1);
    
    % extract segments
    disp('Extracting segments');
    for segIdx = 1:numSegs
        segId = sortedIds(segIdx);
        segBox = segProps(segIdx).BoundingBox;
        
        % show progress
        Util.showProgress( ...
            segIdx, numSegs, ['Segment ', num2str(segId)]);
        
        % ATTENTION!
        % just as isosurface, also regionprops expects the
        % following axis to dimension mapping:
        % Y is dim. one, X is dim. two, and Z is dim. three
        boxMin = round([segBox(2), segBox(1), segBox(3)]);
        boxWidth = [segBox(5), segBox(4), segBox(6)];
        boxMax = boxMin + boxWidth - 1;
        
        % build bounding boxes
        boxLocal = [boxMin(:), boxMax(:)];
        boxGlobal = bsxfun( ...
            @plus, boxLocal, cubeBoxBigMin - 1);
        
        % extracting mask
        mask = seg( ...
            boxMin(1):boxMax(1), ...
            boxMin(2):boxMax(2), ...
            boxMin(3):boxMax(3));
        mask = (mask == segIdx);
        
        % store result
        out.segIds(segIdx) = segId;
        out.boxLocal(segIdx, :, :) = boxLocal;
        out.boxGlobal(segIdx, :, :) = boxGlobal;
        out.masks{segIdx} = mask;
    end
    
    % save result
    disp('Storing results...');
    segMaskFile = [cubeFolder, 'segMasks.mat'];
    save(segMaskFile, '-v7.3', '-struct', 'out');

    % show when done
    Util.showProgress( ...
        numSegs, numSegs, 'Done!');
end

function segIds = getSegIds(cubeParams)
    load(cubeParams.segmentFile, 'segments');
    [segIds] = segments.Id;
    
    % make sure segIds are uint32
    segIds = uint32(segIds);
end
