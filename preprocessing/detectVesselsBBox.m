function regionInfo = detectVesselsBBox(inParam,outParam,thisSliceBbox,regionsToRemove)
% function to perform the bv detection on cluster
% assumes wkw data format!

raw = loadRawData(inParam, thisSliceBbox);
vessels = false(size(raw));
for j=1:size(raw,3)
    vessels(:,:,j) = detectVesselsSingleImage(raw(:,:,j));
end
% Regionprops CC of pixel (small compared to blood vessel)
rp = regionprops(vessels, {'Area' 'Centroid' 'PixelIdxList'});
regionInfo = '';
for j=1:length(rp)
    regionInfo = strcat(regionInfo,['Region ' num2str(j) ': ' num2str(rp(j).Area) ' voxel, centroid: ' num2str(round(rp(j).Centroid([2 1 3]))+ thisSliceBbox(:,1)') '\n']);
end
if exist('regionsToRemove','var')
    vessels(cat(1,rp(regionsToRemove).PixelIdxList)) = 0;
end
vessels = imclose(vessels, ones(1, 1, 7));
maskedRaw = raw;
maskedRaw(vessels) = uint8(121);
saveRawData(outParam, thisSliceBbox(:,1)', maskedRaw);
saveSegDataGlobal(struct( ...
    'root', strrep(outParam.root, '/color/', '/segmentation/')), thisSliceBbox(:, 1)', uint32(vessels));
