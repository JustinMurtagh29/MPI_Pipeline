% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

imConfigs = struct;

imConfigs(1).sizeNm = 2000;
imConfigs(1).pos = 1 + [2977, 7204, 1527];
imConfigs(1).view = 'yz';

imConfigs(2).sizeNm = 2000;
imConfigs(2).pos = 1 + [2711, 7348, 1519];
imConfigs(2).view = 'xz';

outSize = [600, 600];

%% Load config data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

%% Generate images
for curConfig = imConfigs
    curSize = round(curConfig.sizeNm ./ param.raw.voxelSize);
    
    curImSize = curSize(ismember('xyz', curConfig.view));
    curSize(not(ismember('xyz', curConfig.view))) = 1;
    
    curBox = round(curConfig.pos - curSize / 2);
    curBox = curBox(:) + [zeros(3, 1), curSize(:) - 1];
    
    curData = loadRawData(param.raw, curBox);
    curData = reshape(curData, curImSize);
    
    curGrid = { ...
        linspace(1, curImSize(2), outSize(2)), ...
        linspace(1, curImSize(1), outSize(1))};
   [curGrid{:}] = ndgrid(curGrid{:});
   
    curData = uint8(interp2( ...
        double(curData), curGrid{:}, 'linear'));
    figure; imshow(curData);
end
