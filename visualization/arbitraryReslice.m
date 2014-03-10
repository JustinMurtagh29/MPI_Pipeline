function arbitraryReslice(raw, seg, mito, v, parameter, id, problemNumber, colors)

% % Not that important, decide which segments to color in video & make sure
% % start and end segment always have the same color by always using setting
% % start segment as object 1 and end object as object 2
segNew = zeros(size(seg));
segNew(seg == v(1).missions(problemNumber).start.id) = 1;
segNew(seg == v(1).missions(problemNumber).end.id) = 2;
% set all possible ends to and id starting with number 3
idx = 3;
for i=1:length(v(1).missions(problemNumber).possibleEnds)
    % If condition to avoid relabeling the end segment as this is also in the
    % possible ends list
    if v(1).missions(problemNumber).possibleEnds(i).id ~= v(1).missions(problemNumber).end.id
        segNew(seg == v(1).missions(problemNumber).possibleEnds(i).id) = idx;
        idx = idx + 1;
    end
end
% % Use only as many colors (of an otherwise persitent colormap) as there are
% % segments colored in this problem 
colors = colors(1:max(segNew(:)),:);

% for each viewport
for vp = 1:length(v)
    % set/extract some needed parameters
    frameSize = [v(vp).width v(vp).height]; % same as in mission.json
    renderedFrames = [40 40]; % setting from levelcreator (how many frames before and after problem location)
    voxelSize = parameter.settings.scale; % from settings.json (anisotropy)
    voxelSizeAfter = [11.24 11.24 11.24]; % size of voxel after interpolation (in video)
    
    % Find out how much data we need around error center (calculations done
    % on small cube for speed purposes)
    border = sqrt(2) .* max(voxelSizeAfter .* [frameSize/2 max(renderedFrames)]); % border needed in nm
    % sqrt(2) ue to maximal effect if reslice is 45 degress wrt original
    % orientation in plane
    border = ceil(repmat(border,1,3) ./ voxelSize); % divide by voxel size of original data and round up

    % Calculate bounding box from border
    bbox = [v(vp).missions(problemNumber).errorCenter - border; v(vp).missions(problemNumber).errorCenter + border]; %bounding box in dataset coordinates
    bbox = bbox - repmat(parameter.bboxBig(:,1) - 1,1,2)'; % coords of segNew/raw voxel (subtract offset of the section in dataset)
    
    % cut local stack for interpolation & rendering according to bbox
    rawLocal = raw(bbox(1,1):bbox(2,1),bbox(1,2):bbox(2,2),bbox(1,3):bbox(2,3));
    segLocal = segNew(bbox(1,1):bbox(2,1),bbox(1,2):bbox(2,2),bbox(1,3):bbox(2,3));
    mitoLocal = mito(bbox(1,1):bbox(2,1),bbox(1,2):bbox(2,2),bbox(1,3):bbox(2,3));


    %% PROBABLY THE INTERESTING PART
    % extract start position & ONB from mission.json (just for shorter
    % notation afterwards)
    direction = v(vp).missions(problemNumber).start.direction([2 1 3]);
    orth1 = v(vp).missions(problemNumber).start.orth1([2 1 3]);
    orth2 = v(vp).missions(problemNumber).start.orth2([2 1 3]);
    % Define grid for original values (rawLocal & segLocal & mitoLocal),
    % which has size 2* border+1 and is centered on the error center
    [X,Y,Z] = meshgrid(-border(1):border(1),-border(2):border(2),-border(3):border(3)); % grid in dataset coordinates [voxel]
    X = X .* voxelSize(1); % grid in dataset coordinates [nm]
    Y = Y .* voxelSize(2);
    Z = Z .* voxelSize(3);
    % Construct rotation matrix from orthogonal basis (one of which is direction)
    A = [orth1; orth2; direction];
    % Define grid for interpolating values (size according to settings from
    % viewport & levelcreator
    [XI, YI, ZI] = meshgrid(-frameSize(1)/2:frameSize(1)/2, -frameSize(2)/2:frameSize(2)/2, -renderedFrames(1):renderedFrames(2)); %new grid [voxel]
    % Linearize representation of points and scale according to desired
    % voxel size
    coords = bsxfun(@times, [XI(:) YI(:) ZI(:)], voxelSizeAfter); % new grid transform to [nm] including making 3D Matrices X,Y,Z to 1D for rotation
    coords = coords*A; % rotate new grid according to ONB
    XI = reshape(coords(:,1),size(XI)); % Undo linearization
    YI = reshape(coords(:,2),size(YI));
    ZI = reshape(coords(:,3),size(ZI));
    % Do the interpolation from original grid (X,Y,Z) indicating location
    % of (raw/seg/mito)Local values onto rotated grid with different
    % anisotropy (all grids are in nm)
    rawLocalNewDirection = interp3(X, Y, Z, rawLocal, XI, YI, ZI, 'cubic');
    segLocalNewDirection = interp3(X, Y, Z, single(segLocal), XI, YI, ZI, 'nearest');
    mitoLocalNewDirection = interp3(X, Y, Z, mitoLocal, XI, YI, ZI, 'nearest');
    % Make movies of stacks (no more processing done in there, just
    % displaying frame by frame and capturing)
    makeTaskMovieSet(rawLocalNewDirection, mitoLocalNewDirection, ['C:\Users\mberning\Desktop\problemInspector\' id v(vp).name 'mito' num2str(v(vp).missions(problemNumber).missionId, '%.3i') '.avi']);
    makeTaskMovieSetColors(uint16(segLocalNewDirection), rawLocalNewDirection, ['C:\Users\mberning\Desktop\problemInspector\' id v(vp).name 'seg' num2str(v(vp).missions(problemNumber).missionId, '%.3i') '.avi'], colors);

end
    
end
