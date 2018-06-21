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

randDends = struct;

randDends(1).title = 'Non-whole cell';
randDends(1).dendIds = find(~ismember( ...
    conn.denMeta.targetClass, {'Somata', 'WholeCell'}));

randDends(2).title = 'Smooth dendrite';
randDends(2).dendIds = find( ...
    conn.denMeta.targetClass == 'SmoothDendrite');

randDends(3).title = 'Apical dendrite';
randDends(3).dendIds = find( ...
    conn.denMeta.targetClass == 'ApicalDendrite');

randDends(4).title = 'Axon initial segment';
randDends(4).dendIds = find( ...
    conn.denMeta.targetClass == 'AxonInitialSegment');

% Shuffle and select a random subset
rng(0);
for curIdx = 1:numel(randDends)
    curDendIds = randDends(curIdx).dendIds;
    
    curDendIds = curDendIds(randperm(numel(curDendIds)));
    curDendIds = curDendIds(1:min(numel(curDendIds), numDendrites));
    
    randDends(curIdx).dendIds = curDendIds;
end

rng(0);
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

for curClassIdx = 1:numel(randDends)
    curClass = randDends(curClassIdx);
    
    curTitle = curClass.title;
    curTag = strrep(lower(curTitle), ' ', '-');
    
    curDendIds = curClass.dendIds;
    curDigits = ceil(log10(1 + numel(curDendIds)));

    for curIdx = 1:numel(curDendIds)
        curDendId = curDendIds(curIdx);
        curParentId = conn.denMeta.parentId(curDendId);
        curAgglo = dendrites(curParentId);

        curSkel = Superagglos.toSkel(curAgglo, skel);
        curSkel.names{1} = sprintf('%s %d', curTitle, curDendId);

        curNmlName = sprintf( ...
            '%0*d_%s-%d.nml', curDigits, curIdx, curTag, curDendId);
        curSkel.write(fullfile(outputDir, curNmlName));
    end
end
