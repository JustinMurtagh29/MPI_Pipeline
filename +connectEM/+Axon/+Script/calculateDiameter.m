% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a.mat');

%% loading data
conn = load(connFile);

param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

segSizes = Seg.Global.getSegToSizeMap(param);

segCentroids = Seg.Global.getSegToCentroidMap(param);
segCentroids = segCentroids .* param.raw.voxelSize;

segCov = load(fullfile(rootDir, 'globalSegmentPCA.mat'), 'covMat');
segCov = reshape(segCov.covMat, [], 3, 3);

%% calculate diameters
[segIds, coms, diams] = Agglo.calculateDiameter( ...
    segSizes, segCentroids, segCov, conn.axons, 'nhoodThresh', 750);