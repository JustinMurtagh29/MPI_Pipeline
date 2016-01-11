function arbitraryResliceAgglo(p, graph, com, idsAgglo, idQuerried, idsMerger, outputFile)

    comStartSegment = com(idQuerried,:);
    % Take 5 segments nearest for direction
    idsBefore = idsAgglo(1:end-1);
    comSegmentsBefore = com(idsBefore,:);
    distances = sqrt(sum(bsxfun(@times, bsxfun(@minus, comSegmentsBefore, comStartSegment), [11.24 11.24 28]).^2,2));
    [sD, perm] = sort(distances, 'ascend');
    nrSeg = min(5, length(perm));
    display(num2str(sD(1:nrSeg)));
    % Direction and ONB for video
    direction = comStartSegment - mean(comSegmentsBefore(perm(1:nrSeg),:),1);
    direction = direction .* [11.24 11.24 28];
    direction = direction ./ norm(direction);
    [orth1, orth2] = findOrthogonals(direction);
    % Define bounding box
    sizeBBox = [500 500 200];
    bboxToLoad(:,1) = comStartSegment - sizeBBox;
    bboxToLoad(:,2) = comStartSegment + sizeBBox;

    % Load raw data and segmentation around comStartSegment
    rawLocal = single(readKnossosRoi(p.raw.root, p.raw.prefix, bboxToLoad));
    segLocal = readKnossosRoi(p.seg.root, p.seg.prefix, bboxToLoad, 'uint32', '', 'raw');
    rawLocal = permute(rawLocal, [2 1 3]);
    segLocal = permute(segLocal, [2 1 3]);

    % Color in only agglomerated segments for now
    segLocalAgglo = zeros(size(segLocal));
    for i=1:length(idsAgglo)
        segLocalAgglo(segLocal == idsAgglo(i)) = 1;
    end
    clear segLocal;

    % set/extract some needed parameters, use B4B settings for now
    frameSize = [640 480];
    renderedFrames = [20 200];
    voxelSize = [11.24 11.24 28];
    voxelSizeAfter = [11.24 11.24 11.24]; % size of voxel after interpolation (in video)

    % Define grid for original values (rawLocal & segLocal),
    % which has size 2* border+1 and is centered on the error center
    [X,Y,Z] = meshgrid(-sizeBBox(1):sizeBBox(1),-sizeBBox(2):sizeBBox(2),-sizeBBox(3):sizeBBox(3)); % grid in dataset coordinates [voxel]
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
    clear coords;
    % Do the interpolation from original grid (X,Y,Z) indicating location of (raw/seg/mito)Local values onto rotated grid with different
    % anisotropy (all grids are in nm)
    rawLocalNewDirection = interp3(X, Y, Z, single(rawLocal), XI, YI, ZI, 'cubic');
    segLocalNewDirection = interp3(X, Y, Z, single(segLocalAgglo), XI, YI, ZI, 'nearest') > 0;
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

