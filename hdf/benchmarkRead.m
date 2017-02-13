function dataSizeOnDisk = benchmarkRead( readFunction, filename, linespec, data )

% Determine file size
[~,dataSizeOnDisk] = system(['du -hs ' filename ' | awk ''NR==1{print $1}''']);

% Test hdf5 reading
sizeRead = 2.^(0:9);
datasetSize = 512;
nrRandomReads = 10;

rng default;
voxelTroughput = cell(1,length(sizeRead));
for i=1:length(sizeRead)
    t = zeros(1,nrRandomReads);
    thisReadSize = repmat(sizeRead(i), 3, 1);
    for j=1:nrRandomReads
        pos = randi(datasetSize - sizeRead(i) + 1, 3,1);
        tic;
        temp = readFunction(filename, pos, thisReadSize);
        t(j) = toc;
    end
    voxelTroughput{i} = prod(thisReadSize) ./ t;
    display(num2str(i));
end

assert(all(data(:) == temp(:)));

for i=1:length(voxelTroughput)
    meanToPlot(i) = mean(voxelTroughput{i});
    stdToPlot(i) = std(voxelTroughput{i});
end

errorbar(meanToPlot, stdToPlot, linespec);

end

