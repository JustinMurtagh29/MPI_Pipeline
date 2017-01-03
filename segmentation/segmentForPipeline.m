function segmentForPipeline( root, prefix, bbox, segFunction, saveFile )

% Load classification
aff = loadClassData(struct('root', root, 'prefix', prefix), bbox);
aff = imcomplement(aff);

% Perform segmentation
seg = segFunction(aff);
seg = uint16(seg{1,1});

% If folder does not exist, create it
saveFolder = fileparts(saveFile);
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end

% Save segmentation to MATLAB file in 'saveFolder'
Util.save(saveFile, seg);

end

