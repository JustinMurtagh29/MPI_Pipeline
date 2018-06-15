% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');
outputDir = '/home/amotta/Desktop/neurites-for-path-length-calibration';

numAxons = 50;
numDendrites = 50;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

conn = connectEM.Connectome.load(param, connFile);

dendrites = load(conn.info.param.dendriteFile);
dendrites = dendrites.dendrites;

voxelSize = param.raw.voxelSize;
points = Seg.Global.getSegToPointMap(param);

%% Select neurites
rng(0);

% For whole cells we have proper skeleton tracings to do the calibration,
% so we are skipping somata and whole cells here.
randDendIds = find(~ismember( ...
    conn.denMeta.targetClass, {'Somata', 'WholeCell'}));
randDendIds = randDendIds(randperm(numel(randDendIds), numDendrites));

randAxonIds = randperm(numel(conn.axons), numAxons);
randAxonIds = reshape(randAxonIds, [], 1);

%% Generate skeletons
skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

%% Generate axon NMLs
clear cur*;
curDigits = ceil(log10(1 + numAxons));

for curIdx = 1:numAxons
    curAxonId = randAxonIds(curIdx);
    curAgglo = SuperAgglo.fromAgglo( ...
        conn.axons(curAxonId), points, ...
        'mst', 'voxelSize', voxelSize);
    
    curSkel = Superagglos.toSkel(curAgglo, skel);
    curSkel.names{1} = sprintf('Axon %d', curAxonId);
    
    curNmlName = sprintf( ...
        '%0*d_axon-%d.nml', ...
        curDigits, curIdx, curAxonId);
    curSkel.write(fullfile(outputDir, curNmlName));
end

%% Generate dendrites NMLs
clear cur*;
curDigits = ceil(log10(1 + numDendrites));

for curIdx = 1:numDendrites
    curDendriteId = randDendIds(curIdx);
    curParentId = conn.denMeta.parentId(curDendriteId);
    curAgglo = dendrites(curParentId);
    
    curSkel = Superagglos.toSkel(curAgglo, skel);
    curSkel.names{1} = sprintf('Dendrite %d', curDendriteId);
    
    curNmlName = sprintf( ...
        '%0*d_dendrite-%d.nml', ...
        curDigits, curIdx, curDendriteId);
    curSkel.write(fullfile(outputDir, curNmlName));
end
