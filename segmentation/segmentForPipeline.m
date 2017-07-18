function segmentForPipeline(affParam, bbox, segFunction, saveFile)

% Load classification
aff = loadClassData(affParam, bbox);
mask = aff ~= -2;  % get mask for nonpadded areas
aff = imcomplement(aff);

% Perform segmentation
seg = segFunction(aff);
seg = uint16(seg{1,1});

% Remove small objects (less than 9 voxels) from segmentation
segMask = bwareaopen(seg, 9);
clear seg;
affImposed = imimposemin(aff, segMask);
clear segMask;
seg = uint16(watershed(affImposed, 26));
clear affImposed;
seg(~mask) = 0;  %delete segments in the padded area

% If folder does not exist, create it
saveFolder = fileparts(saveFile);
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end

% Save segmentation to MATLAB file in 'saveFolder'
Util.save(saveFile, seg);

end

