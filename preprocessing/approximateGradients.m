function [raw_mean, x, y, z] = approximateGradients(vesselsMasked, bbox, filterSize)

% Make multiple of filter size by removing on upper limit
upperLimitCapped = floor((bbox(:,2)-bbox(:,1)+1)./filterSize) .* filterSize - 1 + bbox(:,1);
bbox(:,2) = upperLimitCapped;

% For generating grid for interpolation in original coordinates
x = (bbox(1,1)+filterSize(1)/2-0.5):filterSize(1):(bbox(1,2)-filterSize(1)/2+0.5);
y = (bbox(2,1)+filterSize(2)/2-0.5):filterSize(2):(bbox(2,2)-filterSize(2)/2+0.5);
z = (bbox(3,1)+filterSize(3)/2-0.5):filterSize(3):(bbox(3,2)-filterSize(3)/2+0.5);

% Mean downsampling
display('Mean downsampling');
tic;
zCoords = bbox(3,1):bbox(3,2);
zCoords = mat2cell(zCoords, 1, repmat(filterSize(3), (bbox(3,2)-bbox(3,1)+1)/filterSize(3), 1));
thisSliceBbox = bbox;
for i=1:length(zCoords)
    thisSliceBbox(3,:) = [zCoords{i}(1) zCoords{i}(end)];
    raw = loadRawData(vesselsMasked, thisSliceBbox);
    raw_mean(:,:,i) = nlfilter3(raw, @mean, filterSize);
    Util.progressBar(i, length(zCoords));
end
toc;

end
