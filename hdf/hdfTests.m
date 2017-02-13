%% Load data
% load /gaba/u/mberning/results/pipeline/20151111T183414/allParameter.mat;
% data = load1024SegCube(p);
% save('/home/mberning/Desktop/segData.mat', 'data', '-v7.3');

%% Load local version instead
load('/home/mberning/Desktop/segData.mat');

%% Settings
testFilename1 = '/home/mberning/Desktop/test_nocomp.hdf';
testFilename2 = '/home/mberning/Desktop/test_comp.hdf';
testFilename3 = '/home/mberning/Desktop/testKhierachy/';
% testFilename4 = '/home/mberning/Desktop/testMhierachy/';

%% Write different encodings of data
h5create(testFilename1, '/seg', size(data), 'Datatype', 'uint32', 'ChunkSize', [32 32 32], 'Deflate', 0);
h5write(testFilename1, '/seg', data);
h5create(testFilename2, '/seg', size(data), 'Datatype', 'uint32', 'ChunkSize', [32 32 32], 'Deflate', 1);
h5write(testFilename2, '/seg', data);
writeKnossosRoi(testFilename3, 'seg', [1 1 1], data, 'uint32');
% mortonSaveRoi(testFilename4, 'seg', [1 1 1], data)

%% Perform benchmark & plot
figure;
hold on;
benchmarkRead(@(x,y,z)h5read(x,'/seg', y, z), testFilename1, '-b', data);
benchmarkRead(@(x,y,z)h5read(x,'/seg', y, z), testFilename2, '-r', data);
benchmarkRead(@(x,y,z)readKnossosRoi(x, 'seg', [y, y + z - 1], 'uint32'), testFilename3, '-g', data);
% benchmarkRead(@(x,y,z)mortonLoadRoi(x, 'seg', [y; z]), testFilename4);

% Make nicer
ylim([0 1e9]);
xlabel('Side of loaded cube [voxel]');
ylabel('Voxel troughput [voxel/s]');
legend('HDF5 high-level matlab; no compression', 'HDF5 high-level matlab; minimal gzip compression (200MB)', 'KNOSSOS hierachy, uncompressed');
set(gca, 'YScale', 'log');
set(gca, 'XTick', 1:11);
set(gca, 'XTickLabel', strseq('', 2.^(0:10)));
