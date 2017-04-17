function score = findQueries(agglo, segmentMeta)
voxelRepresentation = connectEM.growout_mag4(agglo, segmentMeta);
bwconncomp(voxelRepresentation)
ballF1 = @(x,y,z)cat(4, repmat((-x:x)',1,2*y+1,2*z+1)/x, repmat((-y:y),2*x+1,1,2*z+1)/y, repmat(reshape(-z:z,[1,1,2*z+1]),2*x+1, 2*y+1,1)/z);
ballF2 = @(x,y,z)sqrt(sum(ballF1(x,y,z).^2,4)) <= 1;
connectEM.makesurf(voxelRepresentation, '1.issf');
voxelRepresentation = imopen(voxelRepresentation, ballF2(2, 2, 1)); %to get rid of spines necks
xx = 1
connectEM.makesurf(voxelRepresentation,'2.issf');
xx = 2
%conncomps_super = bwconncomp(voxelRepresentation);
%[~, idx_max] = max(cellfun(@length, conncomps_super.PixelIdxList));
%voxelRepresentation(setdiff(1 : numel(voxelRepresentation), conncomps_super.PixelIdxList{idx_max})) = false; % remove spine heads
%connectEM.makesurf(voxelRepresentation, '3.issf');
voxelRepresentation = imclose(voxelRepresentation, ballF2(15, 15, 6)); % to close gaps in agglomeration
xx = 4
connectEM.makesurf(voxelRepresentation, '4.issf');
xx = 5

surround = [50, 50, 20];
voxelRepresentation = padarray(voxelRepresentation, ...
    surround);

score = zeros(size(voxelRepresentation));
for idx1 = 1 + surround(1) : size(voxelRepresentation, 1) - surround(1)
    for idx2 = 1 + surround(2) : size(voxelRepresentation, 2) - surround(2)
        for idx3 = 1 + surround(3) : size(voxelRepresentation, 3) - surround(3)
            if voxelRepresentation(idx1, idx2, idx3)
                localSurround = voxelRepresentation( ...
                    idx1 - surround(1) : idx1 + surround(1), ...
                    idx2 - surround(2) : idx2 + surround(2), ...
                    idx3 - surround(3) : idx3 + surround(3));
                conncomps = bwconncomp(localSurround);
                centervoxelIdx = sub2ind(size(localSurround), surround(1) + 1, surround(2) + 1, surround(3) + 1);
                idx_center = find(cellfun(@(x)ismember(centervoxelIdx, x), conncomps.PixelIdxList));
                [coordsX, coordsY, coordsZ] = ind2sub(conncomps.ImageSize, conncomps.PixelIdxList{idx_center});
                list_coords = [coordsX, coordsY, coordsZ];
                my_pca = pca(list_coords);
                score(idx1, idx2, idx3) = mean(list_coords * my_pca(:, 1));
            end
        end
    end
end
