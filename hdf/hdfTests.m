%% Load data
load /gaba/u/mberning/results/pipeline/20151111T183414/allParameter.mat;
data = load1024SegCube(p);
save('/home/mberning/Desktop/segData.mat', 'data', '-v7.3');

%% Load local version instead
load('/home/mberning/Desktop/segData.mat');

%% Settings
testFilename = '/home/mberning/Desktop/test.hdf';

%% Test hdf5 writing
h5create(testFilename, '/seg', size(data), 'Datatype', 'uint32', 'ChunkSize', [32 32 32], 'Deflate', 0);
h5write(testFilename, '/seg', data);

%% Determine file size
fileInfo = dir(testFilename);

%% Test hdf5 reading
sizeRead = 2.^(0:10);
datasetSize = 1024;
nrRandomReads = 10;

rng default;
voxelTroughput = cell(1,length(sizeRead));
for i=1:length(sizeRead)
    t = zeros(1,nrRandomReads);
    thisReadSize = repmat(sizeRead(i), 3, 1);
    for j=1:nrRandomReads
        pos = randi(datasetSize - sizeRead(i) + 1, 3,1);
        tic;
        data = h5read(testFilename,'/seg', pos, thisReadSize);
        t(j) = toc;
    end
    voxelTroughput{i} = prod(thisReadSize) ./ t;
    display(num2str(i));
end

for i=1:length(voxelTroughput)
    meanToPlot(i) = mean(voxelTroughput{i});
    stdToPlot(i) = std(voxelTroughput{i});
end

errorbar(meanToPlot, stdToPlot);
ylim([0 12e8]);
ylabel('Voxel troughput [voxel/s]');
xlabel('Side of loaded cube [voxel]');
legend('HDF5 high-level matlab');
