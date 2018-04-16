% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
outDir = '/tmpscratch/amotta/l4/2018-04-16-distance-volume-test/wkw';

distRefIsoFile = '/tmpscratch/amotta/l4/2018-04-11-smooth-dendrite-isosurfaces/mat/iso-9.mat';
distThreshUm = 10;

boxSize = [128, 128, 128];

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

wkwInit('new', outDir, 32, 32, 'uint32', 1);

%% Build boxes
boxGlobal = param.bbox;
assert(~any(mod((boxGlobal(:, 1)' - 1), boxSize)));

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

%%
taskArgs = cellfun(@(box) {{box}}, boxes);
taskSharedArgs = {param, outDir, distRefIsoFile};

job = Cluster.startJob( ...
    @taskFunction, taskArgs, ...
    'sharedInputs', taskSharedArgs, ...
    'cluster', cluster, ...
    'name', mfilename);

%% Utility
function taskFunction(param, outDir, distRefIsoFile, bbox)
    import connectEM.Availability.buildDistanceVolume;
    
    distRefIso = load(distRefIsoFile);
    distRefIso = distRefIso.isoSurf;
    
    distVol = buildDistanceVolume(param, distRefIso, bbox);
    distVol = uint32(distVol);
    
    wkwSaveRoi(outDir, bbox(:, 1)', distVol);
end
