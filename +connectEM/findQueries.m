function [voxelRepresentation3, mincoord] = findQueries(agglo, segmentMeta, options)
[voxelRepresentation, mincoord] = connectEM.growout_mag4(agglo, segmentMeta);

ballF1 = @(x,y,z)cat(4, repmat((-x:x)',1,2*y+1,2*z+1)/x, repmat((-y:y),2*x+1,1,2*z+1)/y, repmat(reshape(-z:z,[1,1,2*z+1]),2*x+1, 2*y+1,1)/z);
ballF2 = @(x,y,z)sqrt(sum(ballF1(x,y,z).^2,4)) <= 1;
% remove anisotropy
voxelRepresentation2 = voxelRepresentation(:, :, imresize(1 : size(voxelRepresentation, 3), [1, round(size(voxelRepresentation, 3) * 2.313/2)*2], 'nearest'));
if options.debug
    connectEM.makesurf(nlfilter3(voxelRepresentation2, @mode, [2, 2, 2]), '1.issf');
end
voxelRepresentation2 = imopen(voxelRepresentation2, ballF2(2, 2, 2)); %to get rid of spines necks
if options.debug
    connectEM.makesurf(nlfilter3(voxelRepresentation2, @mode, [2, 2, 2]),'2.issf');
end
conncomps_super = bwconncomp(voxelRepresentation2);
[~, idx_max] = max(cellfun(@length, conncomps_super.PixelIdxList));
for idx = setdiff(1 : length(conncomps_super.PixelIdxList), idx_max)
    voxelRepresentation2(conncomps_super.PixelIdxList{idx}) = false; % remove spine heads
end
if options.debug
    connectEM.makesurf(nlfilter3(voxelRepresentation2, @mode, [2, 2, 2]), '3.issf');
end
%downsample to mag8
voxelRepresentation2 = nlfilter3(voxelRepresentation2, @mode, [2, 2, 2]);

voxelRepresentation3 = imfilter(double(voxelRepresentation2), double(ballF2(6, 6, 6)));
if options.debug
    save('voxelRepresentation3', 'voxelRepresentation3', '-v7.3');
    connectEM.makesurf(voxelRepresentation2, '4.issf');
    save('voxelRepresentation2', 'voxelRepresentation2', '-v7.3');
end