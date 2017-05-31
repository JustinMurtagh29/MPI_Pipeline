display('Gradient correction');
thisSliceBbox = dataset.bbox;
[X,Y,Z] = meshgrid(y,x,z);
for i=2:9
    % Determine bounding box and interpolation for this z-slice
    thisSliceBbox(3,:) = [zCoords{i}(1) zCoords{i}(end)];	
	xq = thisSliceBbox(1,1):thisSliceBbox(1,2);
    yq = thisSliceBbox(2,1):thisSliceBbox(2,2);
    zq = thisSliceBbox(3,1):thisSliceBbox(3,2);
    [Xq,Yq,Zq] = meshgrid(yq,xq,zq);
    % Read original data and detected vessel
    % Interpolate correction voxel for each voxel and multiply
    correctionForSlice =  interp3(X,Y,Z,correctionVolume,Xq,Yq,Zq, 'linear', 121);
    raw = readKnossosRoi(vesselsMasked.root, vesselsMasked.prefix, thisSliceBbox); 
    vessels =  readKnossosRoi(strrep(vesselsMasked.root, '/color/', '/segmentation/'), vesselsMasked.prefix, thisSliceBbox, 'uint32');
    raw = uint8(correctionForSlice .* double(raw));
	raw(vessels > 0) = 128;
    % Save to new datset to be used in pipeline
    writeKnossosRoi(gradientCorrected.root, gradientCorrected.prefix, thisSliceBbox(:,1)', raw);
    clear raw vessels;
    Util.progressBar(i, length(zCoords));
end
clear X Xq Y Yq Z Zq;