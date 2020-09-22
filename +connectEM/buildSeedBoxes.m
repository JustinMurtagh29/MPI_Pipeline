% This script produces boxes for the unbiased seeding of calibration axons
% following the same strategy used in Boergens, Berning et al. (2017) Nat
% Methods. In brief, N boxes of (2.5 µm)³ are randomly placed within a
% cuboid of (15 µm)³ centered on the center of the EM dataset.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
% Modified here by
%   Sahil Loomba <sahil.loomba@brain.mpg.de>

%% Configuration
param = p;
outDir = fullfile(param.saveFolder,'tracings', 'box-seeded');
mkdir(outDir)
numBoxes = 150;
boxDim = 5000;% 2500 for agglo eval
centerBoxDim = 100000; % 15000 for agglo eval
outFiles = arrayfun(@(x) fullfile(outDir,sprintf('calibration-box-spines-%02d.nml',x)),1:numBoxes,'uni',0);
%% Build boxes
centerPos = round(mean(param.bbox, 2));

% for center of dataset for agglo evaluation
%centerBoxSize = round(15000 ./ param.raw.voxelSize(:));
% for spine-head training data boxes
centerBoxSize = round(centerBoxDim ./ param.raw.voxelSize(:));
centerBox = centerPos(:) - round(centerBoxSize(:) / 2);
centerBox = centerBox + [[0; 0; 0], centerBoxSize - 1];

boxSize = round(boxDim ./ param.raw.voxelSize(:));
boxes = nan(3, 2, numBoxes);

rng(0);
for curIdx = 1:numBoxes
    curOff = centerBoxSize - boxSize + 1;
    curOff = arrayfun(@(n) randi([0, n]), curOff);
    
    curBox = curOff + [[0; 0; 0], boxSize - 1];
    curBox = centerBox(:, 1) + curBox;
    boxes(:, :, curIdx) = curBox;

    % export to skel
    skel = skeleton();
    skel = skel.setParams(param.experimentName,param.raw.voxelSize,[0,0,0]);
    skel = skel.addTree(outFiles{curIdx}(1:end-4), reshape(curBox(:,1),1,''));
    skel = setBoundingBox(skel, Util.convertMatlabToWKBbox(curBox),'user');
    skel.write(outFiles{curIdx});
end

%% Export boxes in webKnossos format
clear cur*;

curDigits = ceil(log10(1 + numBoxes));

for curIdx = 1:numBoxes
    curBox = boxes(:, :, curIdx);
    curBox(:, 2) = 1 + diff(curBox, 1, 2)';
    curBox(:, 1) = curBox(:, 1) - 1;
    
    fprintf( ...
        'Box %0*d: %d, %d, %d, %d, %d, %d\n', ...
        curDigits, curIdx, curBox);
end

