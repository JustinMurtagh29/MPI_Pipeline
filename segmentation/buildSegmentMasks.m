function jobSegMasks(cubeParam)
    function segIds = getSegIds(cubeFolder)
        globalSegFile = [cubeFolder, 'segGlobal.mat'];
        
        % load global seg matrix
        fileStruct = load(globalSegFile, 'seg');
        segMat = fileStruct.seg;
        
        % find unique ids
        uniqueIds = unique(segMat(:));
        segIds = sort(uniqueIds);
    end
    
    function saveSegmentIds(segsDir, segList)
        segsFileStruct = struct;
        segsFileStruct.segments = segList;
        
        segsFile = [segsDir, 'segments.mat'];
        save(segsFile, '-struct', 'segsFileStruct');
    end
        
    cubeFolder = cubeParam.saveFolder;
    cubeBoxBig = cubeParam.bboxBig;
    cubeBoxBigMin = cubeBoxBig(:, 1);
    
    % build path to segments folder
    cubeSegsFolder = [cubeFolder, 'segments', filesep];
    
    % create folder if it does not exist
    if ~exist(cubeSegsFolder, 'dir')
        % create segments folder
        mkdir(cubeSegsFolder);
    end
    
    % load global seg matrix
    disp('Loading global segment matrix');
    sortedIds = getSegIds(cubeFolder);
    
    % save list of segments
    saveSegmentIds(cubeSegsFolder, sortedIds);
    
    % find extremal segment ids
    minId = sortedIds(1);
    maxId = sortedIds(end);
    numSegs = maxId - minId + 1;
    
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
        boxLocal = [boxMin', boxMax'];
        boxGlobal = boxLocal + repmat(cubeBoxBigMin, [1, 2]) - 1;
        
        % extracting mask
        mask = seg( ...
            boxMin(1):boxMax(1), ...
            boxMin(2):boxMax(2), ...
            boxMin(3):boxMax(3));
        mask = (mask == segIdx);
        
        % prepare results
        maskStruct = Util.buildMaskStruct( ...
            boxLocal, boxGlobal, mask);
        
        % build path to results
        segDir = [cubeSegsFolder, num2str(segId), filesep];
        
        % create folder if needed
        if ~exist(segDir, 'dir')
            mkdir(segDir);
        end
        
        maskFile = [segDir, 'mask.mat'];
        save(maskFile, '-struct', 'maskStruct');
    end

    % show when done
    Util.showProgress( ...
        numSegs, numSegs, 'Done!');
end
