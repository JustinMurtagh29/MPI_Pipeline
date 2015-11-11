% Where to find images
pathToImages = '/run/media/mberning/Seagate Expansion Drive/8-bit aligned TIFF/';
% Where to find function readKnossosRoi/Cube functions
addpath(genpath('/home/mberning/fsHest/Personal/berningm/backup/20150120Laptop/code/auxiliaryMethods/'));
% Find all files (will be placed in cubes in alphabetical order)
files = dir([pathToImages '*.tif']);
% Define distribution images are to be fitted to
wantedDistribution = round(10000*normpdf(-127:128, 0, 50));
% Define how many images (z-slices) to load at once
slices = [0:64:length(files) length(files)];
% Define which resolutions to write
resolutions = 2.^(0:6);

%% Do for all files (read images, write KNOSSOS hierachies)
for s=1:length(slices)-1
    % Display progress to command line
    display(['Progress: ' num2str(s) '/' num2str(length(slices)-1)]);
    tic;
    % Preallocate array of 128 z slices as RAM based temp storage
    cube = zeros(9706,10765,slices(s+1)-slices(s),'uint8');
    % Read images
    for i=1:slices(s+1)-slices(s)
        img = imread([pathToImages files(slices(s)+i).name]);
        cube(:,:,i) = histeq(img, wantedDistribution);
    end
    clear img;
    % Write all different Knossos hierachies (different subsamplings)
    for r=1:length(resolutions)
        % Where to store different resolutions
        directory = ['/home/mberning/localStorage/data/stackED/' num2str(resolutions(r)) '/'];
        if ~exist(directory, 'dir'); mkdir(directory); end;
        % Write data
        writeKnossosRoi(directory, ['No14-020515_mag' num2str(resolutions(r))], [1 1 1+(slices(s)/resolutions(r))], cube);
        % Subsample for next resolution
        if r ~= length(resolutions)
            sC = size(cube);
            cube = cube(1:end-mod(sC(1),2),1:end-mod(sC(2),2),1:end-mod(sC(3),2));
            cube = nlfilter3(cube, @mean, [2 2 2]);
        end
    end
    toc;
end
