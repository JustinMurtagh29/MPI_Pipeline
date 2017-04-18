function score = findQueries(agglo, segmentMeta)
voxelRepresentation = connectEM.growout_mag4(agglo, segmentMeta);
bwconncomp(voxelRepresentation)
ballF1 = @(x,y,z)cat(4, repmat((-x:x)',1,2*y+1,2*z+1)/x, repmat((-y:y),2*x+1,1,2*z+1)/y, repmat(reshape(-z:z,[1,1,2*z+1]),2*x+1, 2*y+1,1)/z);
ballF2 = @(x,y,z)sqrt(sum(ballF1(x,y,z).^2,4)) <= 1;
voxelRepresentation2 = voxelRepresentation(:, :, imresize(1 : size(voxelRepresentation, 3), [1, size(voxelRepresentation, 3) * 2.313], 'nearest'));
%connectEM.makesurf(voxelRepresentation2, '1.issf');
tic
voxelRepresentation2 = imopen(voxelRepresentation2, ballF2(2, 2, 2)); %to get rid of spines necks
toc
xx = 1
tic
%connectEM.makesurf(voxelRepresentation2,'2.issf');
toc
xx = 2
tic
conncomps_super = bwconncomp(voxelRepresentation2);
toc
[~, idx_max] = max(cellfun(@length, conncomps_super.PixelIdxList));
tic
for idx = setdiff(1 : length(conncomps_super.PixelIdxList), idx_max)
    voxelRepresentation2(conncomps_super.PixelIdxList{idx}) = false; % remove spine heads
end
toc

%connectEM.makesurf(voxelRepresentation2, '3.issf');
tic
voxelRepresentation2 = imclose(voxelRepresentation2, ballF2(4, 4, 4)); % to close gaps in agglomeration
toc
xx = 4
connectEM.makesurf(voxelRepresentation2, '4.issf');
xx = 5
save('voxelRepresentation2', 'voxelRepresentation2', '-v7.3');
xx = 6
voxelRepresentation2 = voxelRepresentation2 > 0.1;
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
