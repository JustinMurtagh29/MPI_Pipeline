function writeSupercubes(root, prefix, bbox, magToWrite, outputDirectory)
% Write subsampled data, it is assumed root points to magnification one power of 2 less
% than the first mag to write & all magnifications increase by a power of 2

% Read data once
raw = readKnossosRoi(root, prefix, bbox);

% Write multiple resoultions
for i=1:length(magToWrite)
    if size(raw,3) > 1
        raw = nlfilter3(raw, @median, [2 2 2]);
        [startIdx, endIdx] = regexp(prefix, '_mag.{1,3}');
        newPrefix = [prefix(1:startIdx-1) '_mag' num2str(magToWrite(i))];
        startPos = (bbox(:,1)' - [1 1 1]) ./ magToWrite(i) + [1 1 1];
        writeKnossosRoi(fullfile(outputDirectory, num2str(magToWrite(i))), newPrefix, startPos, raw, 'uint8', '', 'noRead');
    end
end

end

