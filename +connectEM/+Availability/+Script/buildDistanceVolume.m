% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
outDir = '/tmpscratch/amotta/l4/2018-04-16-distance-volume-test/wkw';

distRefIsoFile = '/tmpscratch/amotta/l4/2018-04-11-smooth-dendrite-isosurfaces/mat/iso-9.mat';
distThreshUm = 10;

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

wkwInit('new', outDir, 32, 32, 'uint8', 1);

%% Starting job
cluster = Cluster.getCluster( ...
    '-p 0', ...
    '-l h_vmem=12G', ...
    '-l h_rt=6:00:00');

taskArgs = arrayfun(@(i) {{i}}, 1:numel(param.local));
taskSharedArgs = {param, outDir, distRefIsoFile, distThreshUm};

job = Cluster.startJob( ...
    @taskFunction, taskArgs, ...
    'sharedInputs', taskSharedArgs, ...
    'cluster', cluster, ...
    'name', mfilename);

%% Utility
function taskFunction(param, outDir, distRefIsoFile, distThreshUm, cubeIdx)
    import connectEM.Availability.buildDistanceVolume;
    
    distRefIso = load(distRefIsoFile);
    distRefIso = distRefIso.isoSurf;
    
    bbox = param.local(cubeIdx).bbox;
    
    mask = buildDistanceVolume(param, distRefIso, distThreshUm, bbox);
    if ~any(mask(:)); return; end
    
    mask = uint8(mask);
    wkwSaveRoi(outDir, bbox(:, 1)', mask);
end