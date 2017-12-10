% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
outputDir = '/home/amotta/Desktop/triplets';

chiasmataFile = fullfile( ...
    rootDir, 'tripletDetection', ...
    '20171209T164223-on-axons-13a', ...
    '20171209T164745_chiasmata.mat');

exportAxonIds = [143; 20701; 6426];

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

chiasmata = load(chiasmataFile);
axonFile = chiasmata.info.param.axonFile;
chiasmata = chiasmata.chiasmata;

axons = load(axonFile, 'axons', 'indBigAxons');
axonIds = find(axons.indBigAxons);
axons = axons.axons(axonIds);

% translate axon IDs
[~, exportAxonIdsRel] = ismember(exportAxonIds, axonIds);
assert(all(exportAxonIdsRel));

%% export skeletons
mkdir(outputDir);

for curIdx = 1:numel(exportAxonIds)
    curAxonId = exportAxonIds(curIdx);
    curAxonIdRel = exportAxonIdsRel(curIdx);
    
    curSkel = ...
        connectEM.Chiasma.Detect.buildSkeleton( ...
            axons(curAxonIdRel), chiasmata{curAxonIdRel});
	curSkel = Skeleton.setParams4Pipeline(curSkel, param);
    
	curSkelName = sprintf('%d_axon-%d.nml', curIdx, curAxonId);
    curSkel.write(fullfile(outputDir, curSkelName));
end