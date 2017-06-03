function newPos = correctQueryLocationToEndOfSegment(p, cc, pos, dir, nmStepback)
% Correct query location to center of mass of last ~1000nm of cc along direction vector

    % Define area to load
    sizeBBox = [100 100 40];
    scale = p.raw.voxelSize;
    bboxToLoad(:,1) = pos - sizeBBox;
    bboxToLoad(:,2) = pos + sizeBBox;

    % Load data and set logical indicating whether voxel is in cc
    seg = readKnossosRoi(p.seg.root, p.seg.prefix, bboxToLoad, 'uint32', '', 'raw');
    segAgglo = imdilate(ismember(seg, cc), ones(3,3,3));
    clear seg;

    % Get all true voxel postions, offset to global coordinates
    [pixelList(:,1), pixelList(:,2), pixelList(:,3)] = ind2sub(size(segAgglo), find(segAgglo));
    pixelList = bsxfun(@plus, pixelList, bboxToLoad(:,1)' - [1 1 1]);
    clear props;

    % Project onto direction vector and keep only voxel in last 10 frames
    dirNm = (dir .* scale)' ./ norm(dir .* scale);
    prOnDir = bsxfun(@times, pixelList, scale) * dirNm;
    lastVoxelPosAlongDir = max(prOnDir);
    withinStartRegion = prOnDir > (lastVoxelPosAlongDir - nmStepback - 50) ...
                      & prOnDir < (lastVoxelPosAlongDir - nmStepback + 50);

    % New position is center of mass of those pixel
    newPos = round(mean(pixelList(withinStartRegion,:)));

    % Check whether start is in right process
    %testPos = newPos - bboxToLoad(:,1)' + [1 1 1];
    %if ~segAgglo(testPos(1),testPos(2),testPos(3))
    %    testPos
    %    error('Start position not in agglomeration');
    %end

end

