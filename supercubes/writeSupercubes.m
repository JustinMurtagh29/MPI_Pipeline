function writeSupercubes(root, prefix, bbox, magToWrite, outputDirectory)
% Write subsampled data, it is assumed root points to magnification one power of 2 less
% than the first mag to write & all magnifications increase by a power of 2

% Read data once
raw = readKnossosRoi(root, prefix, bbox);

% Write multiple resoultions
for i=1:length(magToWrite)
    raw = imresize(raw, 0.5);
    newPrefix = [prefix(1:endIdx-1) '_mag' num2str(i)];
    writeKnossosRoi(fullfile(outputDirectory, num2str(magToWrite)), newPrefix, double(bbox(:,1)'), raw);
end

end

