function arbitraryResliceAgglo(p, graph, com, idsAgglo, idQuerried, idsMerger, outputFile)

    comStartSegment = com(idQuerried,:);
    % Take 5 segments before for direction (less if start of agglo)
    nrPreviousSegments = min(length(idsAgglo)-1, 6);
    comOf5SegmentsBefore = com(idsAgglo(end-nrPreviousSegments:end-1),:);
    direction = comStartSegment - comOf5SegmentsBefore;
    direction = direction ./ norm(direction);
    [orth1, orth2] = findOrthogonals(direction);
    sizeBBox = [500 500 200];
    bboxToLoad(:,1) = comStartSegment - sizeBBox;
    bboxToLoad(:,2) = comStartSegment + sizeBBox;

    % Load raw data and segmentation around comStartSegment
    rawLocal = readKnossosRoi(p.raw.root, p.raw.prefix, bboxToLoad);
    segLocal = readKnossosRoi(p.seg.root, p.seg.prefix, bboxToLoad, 'uint32', '', 'raw');

    % Color in only agglomerated segments for now
    segLocalAgglo = zeros(size(seg));
    for i=1:length(idsAgglo)
        segLocalAgglo(segLocal == idsAgglo(i)) = 1;
    end

    % set/extract some needed parameters, use B4B settings for now
    frameSize = [p.kdb.v(2).width p.kdb.v(2).height]; % same as in mission.json
    renderedFrames = [0 200]; % setting from levelcreator (how many frames before and after problem location)
    voxelSize = parameter.settings.scale; % from settings.json (anisotropy)
    voxelSizeAfter = [11.24 11.24 11.24]; % size of voxel after interpolation (in video)

    % Define grid for original values (rawLocal & segLocal),
    % which has size 2* border+1 and is centered on the error center
    [X,Y,Z] = meshgrid(-border(1):border(1),-border(2):border(2),-border(3):border(3)); % grid in dataset coordinates [voxel]
    X = X .* voxelSize(1); % grid in dataset coordinates [nm]
    Y = Y .* voxelSize(2);
    Z = Z .* voxelSize(3);
    % Construct rotation matrix from orthogonal basis (one of which is direction)
    A = [orth1; orth2; direction];
    % Define grid for interpolating values (size according to settings from viewport & levelcreator
    [XI, YI, ZI] = meshgrid(-frameSize(1)/2:frameSize(1)/2, -frameSize(2)/2:frameSize(2)/2, -renderedFrames(1):renderedFrames(2)); %new grid [voxel]
    % Linearize representation of points and scale according to desired voxel size
    coords = bsxfun(@times, [XI(:) YI(:) ZI(:)], voxelSizeAfter); % new grid transform to [nm] including making 3D Matrices X,Y,Z to 1D for rotation
    coords = coords*A; % rotate new grid according to ONB
    XI = reshape(coords(:,1),size(XI)); % Undo linearization
    YI = reshape(coords(:,2),size(YI));
    ZI = reshape(coords(:,3),size(ZI));
    % Do the interpolation from original grid (X,Y,Z) indicating location of (raw/seg/mito)Local values onto rotated grid with different
    % anisotropy (all grids are in nm)
    rawLocalNewDirection = interp3(X, Y, Z, rawLocal, XI, YI, ZI, 'cubic');
    segLocalNewDirection = interp3(X, Y, Z, single(segLocal), XI, YI, ZI, 'nearest');
    % Make movies of stacks (no more processing done in there, just displaying frame by frame and capturing)
    makeTaskMovieSet(rawLocalNewDirection, segLocalNewDirection, outputFile);

end

function [orth1, orth2] = findOrthogonals(v)
    % returns random orthogonal vector and and the 2nd and 3rd vector for ONB 
    v = v ./ norm(v);
    if all(v == [0 0 1])
        orth1 = [1 0 -1].*v([3 2 1]);  
    else
        orth1 = [1 -1 0].*v([2 1 3]);
    end
    orth1 = orth1 ./ norm(orth1);
    orth2 = cross(v, orth1); 
    mix = 2 * (rand(2,1) -0.5);
end

