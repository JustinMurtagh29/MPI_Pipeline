% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
outputDir = '/home/amotta/Desktop/triplets';

chiasmataFile = fullfile( ...
    rootDir, 'tripletDetection', ...
    '20171109T133421-on-axons-10a', ...
    '20171109T133849_chiasmata.mat');

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

chiasmata = load(chiasmataFile);
axonFile = chiasmata.info.param.axonFile;
chiasmata = chiasmata.chiasmata;

axons = load(axonFile, 'axons', 'indBigAxons');
axons = axons.axons(axons.indBigAxons);

%% build chiasmata table
chiasmaT = connectEM.Chiasma.Detect.buildTable(chiasmata, axons);

% look at unsolved triplets
chiasmaT(chiasmaT.isSolved, :) = [];
chiasmaT(chiasmaT.nrExits ~= 3, :) = [];

%% export random agglomerates
rng(0);
randIds = randperm(size(chiasmaT, 1));

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

for curIdx = randIds(1:20)
    curAggloId = chiasmaT.aggloId(curIdx);
    curChiasmaId = chiasmaT.chiasmaId(curIdx);
    
    curSkel = connectEM.Chiasma.Detect.buildSkeleton( ...
        axons(curAggloId), chiasmata{curAggloId}, curChiasmaId);
    curSkel = Skeleton.setParams4Pipeline(curSkel, param);
    
    curFileName = sprintf( ...
        '%02d_axon-%d_chiasma-%d.nml', ...
        curIdx, curAggloId, curChiasmaId);
    curSkel.write(fullfile(outputDir, curFileName));
end