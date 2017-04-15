function score = findQueries(agglo)
voxelRepresentation = growout_mag4(agglo);
ball1 = strel('ball', 5, 2);
ball2 = strel('ball', 15, 6);

voxelRepresentation = imopen(voxelRepresentation, ball1); %to get rid of spines necks
conncomps_super = bwconncomp(voxelRepresentation);
[~, idx_max] = cellfun(@length, conncomps_super.PixelIdxList);
voxelRepresentation(setdiff(1 : numel(voxelRepresentation), conncomps_super.PixelIdxList{idx_max})) = false; % remove spine heads
voxelRepresentation = imclose(voxelRepresentation, ball2); % to close gaps in agglomeration


surround = [50, 50, 20];
voxelRepresentation = paddarray(voxelRepresentation, ...
    surround);

assert(surround - padding <= 0);

score = zeros(size(voxelRepresentation));
for idx1 = 1 + surround(1) : size(voxelRepresentation, 1) - surround(1)
    for idx2 = 1 + surround(2) : size(voxelRepresentation, 2) - surround(2)
        for idx3 = 1 + surround(3) : size(voxelRepresentation, 3) - surround(3)
            if voxelRepresentation(idx1, idx2, idx3)
                localSurround = voxelRepresenstation( ...
                    idx1 - surround(1) : idx1 + surround(1), ...
                    idx2 - surround(2) : idx2 + surround(2), ...
                    idx3 - surround(3) : idx3 + surround(3));
                conncomps = bwconncomp(localSurround);
                centervoxelIdx = sub2ind(size(localSurround), surround(1) + 1, surround(2) + 1, surround(3) + 1);
                [~, idx_center] = cellfun(@(x)ismember(centervoxelIdx, x), conncomps.PixelIdxList);
                [coordsX, coordsY, coordsZ] = ind2sub(conncomps.ImageSize, conncomps.PixelIdxList{idx_center});
                list_coords = [coordsX', coordsY', coordsZ'];
                my_pca = pca(listcoords);
                score(idx1, idx2, idx3) = mean(list_coords * my_pca(:, 1));
            end
        end
    end
end
