function seg20141017( root, prefix, bbox, segFunction, saveFile )

% Load classification
aff = loadClassData(struct('root', root, 'prefix', prefix), bbox);
aff = imcomplement(aff);

% Perform segmentation
seg = segFunction(aff);
seg = uint16(seg{1,1});

% check for mask and, if applicable, apply mask
if exist(fullfile(root,[prefix,'_mask']),'file')
   mask = readKnossosRoi(root, [prefix,'_mask'], bbox, 'logical', '', 'raw');
   seg(~mask) = 0;
end

% If folder does not exist, create it
saveFolder = fileparts(saveFile);
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end

% Save segmentation to MATLAB file in 'saveFolder'
save(saveFile, 'seg');

end

