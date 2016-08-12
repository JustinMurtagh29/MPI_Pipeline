function writeSupercubes(root, prefix, bbox, magToWrite, outputDirectory, subsamplingSegmentation)
    % Write subsampled data, it is assumed root points to magnification one power of 2 less
    % than the first mag to write & all magnifications increase by a power of 2

    % Decide which function to use for downsampling
    % .. and read data once
    if subsamplingSegmentation
        downsamplingFunction = @mode; 
        raw = readKnossosRoi(root, prefix, bbox, 'uint32', '', 'raw');
    else
        downsamplingFunction = @median; 
        raw = readKnossosRoi(root, prefix, bbox);
    end
 
    % Write multiple resoultions
    startPos = bbox(:,1);
    for i=1:length(magToWrite)
        % Make the size of raw divisibale by 2 in each dimension (by padding at upper limit)
        raw(end+1:end+mod(size(raw,1),2),end+1:end+mod(size(raw,2),2),end+1:end+mod(size(raw,3),2)) = 0;
        raw = nlfilter3(raw, downsamplingFunction, [2 2 2]);
        % Deduce name of new file prefix from old file prefix
        [startIdx, endIdx] = regexp(prefix, '_mag.{1,3}');
        newPrefix = [prefix(1:startIdx-1) '_mag' num2str(magToWrite(i))];
        % Update lower coordinate to start writing in same position according to K/wK convention (but different magnification)
        startPos = (startPos - 1) ./ 2 + 1;
        % Write downsampled version
        if subsamplingSegmentation
            writeKnossosRoi(fullfile(outputDirectory, num2str(magToWrite(i))), newPrefix, startPos', raw, 'uint32', '', 'noRead');
        else
            writeKnossosRoi(fullfile(outputDirectory, num2str(magToWrite(i))), newPrefix, startPos', raw, 'uint8', '', 'noRead');
        end
    end

end

