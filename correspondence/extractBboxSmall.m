function seg = extractBboxSmall( segFile, bboxSmall, bboxBig, direction, planesToKeep )
    
    % Construct x, y & z coordinates (voxel global)
    coords{1} = bboxBig(1,1):bboxBig(1,2);
    coords{2}= bboxBig(2,1):bboxBig(2,2);
    coords{3} = bboxBig(3,1):bboxBig(3,2);
    % Decide which part is in bboxSmall (local coordinates of tile)
    idxKeep{1} = ismember(coords{1}, bboxSmall(1,1):bboxSmall(1,2));
    idxKeep{2} = ismember(coords{2}, bboxSmall(2,1):bboxSmall(2,2));
    idxKeep{3} = ismember(coords{3}, bboxSmall(3,1):bboxSmall(3,2));
    if nargin > 3
        % For direction orthogonal to touch of tiles, keep only indicated
        % planes & overwrite earlier indices
        idxKeep{direction} = ismember(coords{direction},planesToKeep);
    end
    load(segFile);
    seg = seg(idxKeep{1},idxKeep{2},idxKeep{3});

end

