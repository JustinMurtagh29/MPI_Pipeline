function gradientCorrectionBBox(inParam,outParam,thisSliceBbox)
% function to perform the gradient correction on cluster
% assumes wkw data format!
datasetVesselsMaskedSeg = inParam;
datasetVesselsMaskedSeg.root = strrep(datasetVesselsMaskedSeg.root, '/color/', '/segmentation/');
vessels =  logical(loadSegDataGlobal(datasetVesselsMaskedSeg, thisSliceBbox));
raw = loadRawData(inParam, thisSliceBbox);

xq = thisSliceBbox(1,1):thisSliceBbox(1,2);
yq = thisSliceBbox(2,1):thisSliceBbox(2,2);
zq = thisSliceBbox(3,1):thisSliceBbox(3,2);
[Xq,Yq,Zq] = meshgrid(yq,xq,zq);
% Read original data and detected vessel
% Interpolate correction voxel for each voxel and multiply
correctionForSlice =  interp3(X,Y,Z,correctionVolume,Xq,Yq,Zq, 'linear', 121);

raw = uint8(correctionForSlice .* double(raw));
raw(vessels > 0) = 121;
% Save to new datset to be used in pipeline
saveRawData(outParam, thisSliceBbox(:, 1)', raw);
