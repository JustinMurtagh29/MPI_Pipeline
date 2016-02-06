function seg20141017( root, prefix, bbox, segFunction, saveFile )

% Load classification
aff = loadClassData(root, prefix, bbox);
aff = imcomplement(aff);

% Perform segmentation
seg = segFunction(aff);
seg = uint16(seg{1,1});

% Save segmentation to MATLAB file in 'saveFolder'
save(saveFile, 'seg');

end

