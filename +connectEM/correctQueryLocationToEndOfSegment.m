function newPos = correctQueryLocationToEndOfSegment(p, cc, pos, dir, nmStepback, visualize)
% Correct query location to point with largest distance to segment surface
% within some bounding box, NOTE: dir and nmStepback not needed anymore, kept
% for compability

    distance = 500;

    % Visualiyation for debug purposes
    if nargin < 6
        visualize = false;
    end

    % Define area to load around border in both directions in each dimension
    surroundToLoad = ceil(repmat(distance,1,3)./p.raw.voxelSize);
    centerOfBbox = surroundToLoad + 1;
    % Bounding box to load from KNOSSOS hierachy centered on position that
    % we want to correct (usually border center of mass right now)
    bboxToLoad(:,1) = pos - surroundToLoad;
    bboxToLoad(:,2) = pos + surroundToLoad;

    % Load segmentation data
    seg = readKnossosRoi(p.seg.root, p.seg.prefix, bboxToLoad, ...
        'uint32', '', 'raw');
    % Calculate distance to watershed border of each voxel in loaded box
    segAggloDistToSurface = bwdist(~ismember(seg, cc));

    % Find all voxel in bbox which would provide maximal node evidence
    maxDistToSurface = sort(unique(segAggloDistToSurface), 'descend');
    display(num2str(maxDistToSurface(1)));
    voxelIds = find(segAggloDistToSurface(:) > (maxDistToSurface(1)-0.1));
    [voxelPos(:,1), voxelPos(:,2), voxelPos(:,3)] = ...
        ind2sub(size(segAggloDistToSurface), voxelIds);

    % Find the one closest to original position if multiple
    distances = pdist2(voxelPos, centerOfBbox);
    [~,idx] = min(distances);

    % New start position in local coordinates
    newPosLocal = voxelPos(idx,:);

    if visualize
        close all;
        figure;
        % Dark grey for agglo segments
        temp = (segAggloDistToSurface > 0) * 0.5;
        % Brighter gray for old start position
        temp(centerOfBbox(1),centerOfBbox(2),centerOfBbox(3)) = 0.75;
        % White for new start position
        temp(newPosLocal(1),newPosLocal(2),newPosLocal(3)) = 1;
        % Play video and pause execution
        implay(temp);
    end

    % New position scale to global offset, also transfer from Matlab
    % coordinates as used here to webKnossos coordinates for API call
    newPos = newPosLocal + bboxToLoad(:,1)' - [2 2 2];

end
