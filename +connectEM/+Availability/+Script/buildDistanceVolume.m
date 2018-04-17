% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
outDir = '/tmpscratch/amotta/l4/2018-04-16-distance-volume-test/wkw';

% Use axon #23325 as reference
% This is the AD-specific axon from panel 4c
distRefIsoFile = '/tmpscratch/amotta/l4/2018-01-24-axons-18a-isosurfaces/mat/iso-23325.mat';

boxSize = [256, 256, 256];

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

distRefPoints = load(distRefIsoFile);
distRefPoints = reducepatch(distRefPoints.isoSurf, 0.1);
distRefPoints = distRefPoints.vertices;

%% Build boxes
boxGlobal = param.bbox;
boxSizeGlobal = 1 + diff(boxGlobal, 1, 2)';
boxCounts = ceil(boxSizeGlobal ./ boxSize);

boxes = cell(prod(boxCounts), 1);
for curIdx = 1:numel(boxes)
   [curX, curY, curZ] = ind2sub(boxCounts, curIdx);
    curBox = ([curX, curY, curZ] - 1) .* boxSize;
    curBox = ([curBox; (curBox + boxSize - 1)])';
    curBox = curBox + boxGlobal(:, 1);
    boxes{curIdx} = curBox;
end

%% Starting job
cluster = Cluster.getCluster( ...
    '-p 0', ...
    '-l h_vmem=12G', ...
    '-l h_rt=6:00:00');

%% Run job
taskArgs = cellfun(@(box) {{box}}, boxes);
taskSharedArgs = {param, outDir, distRefPoints};

% Prepare WKW dataset
wkwInit('new', outDir, 32, 32, 'uint32', 1);

job = Cluster.startJob( ...
    @taskFunction, taskArgs, ...
    'sharedInputs', taskSharedArgs, ...
    'cluster', cluster, ...
    'name', mfilename);

%% Utility
function taskFunction(param, outDir, distRefPoints, bbox)
    import connectEM.Availability.buildDistanceVolume;
    
    distVol = buildDistanceVolume(param, distRefPoints, bbox);
    distVol = uint32(distVol);
    
    wkwSaveRoi(outDir, bbox(:, 1)', distVol);
end
