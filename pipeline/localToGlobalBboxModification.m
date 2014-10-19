function coords = localBboxModification(p)
    % Local segmentation bboxes to use for slicing off overlap for all cubes in volume

    % construct coord array for bBox border area for each cube in local field of p
    coords{1} = [-p.tileBorder(:,1)+[1;1;1] -p.tileBorder(:,1)+[1;1;1]+p.parameter.tileSize];
    coords = repmat(coords,[size(p.local,1),size(p.local,2),size(p.local,3)]);
    % Usually just use inner bounding box, at outer limits of ROI borders needed to always have enough segmentation context
    coords(1,:,:) = cellfun(@(x)([x(1,1)+p.tileBorder(1,1) x(1,2); x(2,1) x(2,2); x(3,1) x(3,2)]),coords(1,:,:), 'UniformOutput', false);
    coords(end,:,:) = cellfun(@(x)([x(1,1) x(1,2)+p.tileBorder(1,2); x(2,1) x(2,2); x(3,1) x(3,2)]),coords(end,:,:), 'UniformOutput', false);
    coords(:,1,:) = cellfun(@(x)([x(1,1) x(1,2); x(2,1)+p.tileBorder(2,1) x(2,2); x(3,1) x(3,2)]),coords(:,1,:), 'UniformOutput', false);
    coords(:,end,:) = cellfun(@(x)([x(1,1) x(1,2); x(2,1) x(2,2)+p.tileBorder(2,2); x(3,1) x(3,2)]),coords(:,end,:), 'UniformOutput', false);
    coords(:,:,1) = cellfun(@(x)([x(1,1) x(1,2); x(2,1) x(2,2); x(3,1)+p.tileBorder(3,1) x(3,2)]),coords(:,:,1), 'UniformOutput', false);
    coords(:,:,end) = cellfun(@(x)([x(1,1) x(1,2); x(2,1) x(2,2); x(3,1) x(3,2)+p.tileBorder(3,2)]),coords(:,:,end), 'UniformOutput', false);


end
