function segmentForPipeline(classParam, bbox, segFunction, saveFile)
    % Load membrane detection
    class = loadClassData(classParam, bbox);

    % Perform segmentation
    seg = segFunction(class);
    seg = uint16(seg);

    % If folder does not exist, create it
    saveFolder = fileparts(saveFile);
    if ~exist(saveFolder, 'dir')
        mkdir(saveFolder);
    end

    % Save segmentation to MATLAB file in 'saveFolder'
    Util.save(saveFile, seg);
end


