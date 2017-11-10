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
axonE3aFile = fullfile( ...
    '/mnt/mpibr/data/Personal/mottaa/L4', ...
    '2017-11-09-Finding-Triplets/data/axons_10_a_plus_allE3a.mat');

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

chiasmata = load(chiasmataFile);
axonFile = chiasmata.info.param.axonFile;
chiasmata = chiasmata.chiasmata;

axons = load(axonFile, 'axons', 'indBigAxons');
axonIds = find(axons.indBigAxons);
axons = axons.axons(axonIds);

% load axons with all E3a flights
axonsE3a = load(axonE3aFile);
axonsE3aChildIds = axonsE3a.childIds;
axonsE3a = axonsE3a.axons;

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

for curIdx = 1:20
    curRandId = randIds(curIdx);
    curAggloId = chiasmaT.aggloId(curRandId);
    curChiasmaId = chiasmaT.chiasmaId(curRandId);
    
    curChiasmata = chiasmata{curAggloId};
    curCenterPos = curChiasmata.nodes( ...
        curChiasmata.ccCenterIdx(curChiasmaId), :);
    
    % axon as before E3a
    curAxon = axons(curAggloId);
    curSkel = connectEM.Chiasma.Detect.buildSkeleton( ...
        curAxon, chiasmata{curAggloId}, curChiasmaId);
    
    % axon as after E3a
    curAxonE3a = axonsE3a(axonsE3aChildIds(axonIds(curAggloId)));
    
    if size(curAxonE3a.nodes, 1) > 1E5
        % restrict to sphere
       [~, curCenterNodeId] = min(pdist2( ...
           curAxonE3a.nodes(:, 1:3) .* param.raw.voxelSize, ...
           curCenterPos .* param.raw.voxelSize, 'squaredeuclidean'));
        curAxonE3a = ...
            connectEM.Chiasma.restrictToSphere( ...
                param, curAxonE3a, curCenterNodeId, 20000);
    end
    
    curSkel = curSkel.addTree( ...
        'After applying all E3a flights', ...
        curAxonE3a.nodes(:, 1:3), curAxonE3a.edges);
    curSkel = Skeleton.setParams4Pipeline(curSkel, param);
    
    curFileName = sprintf( ...
        '%02d_axon-%d_chiasma-%d.nml', ...
        curIdx, curAggloId, curChiasmaId);
    curSkel.write(fullfile(outputDir, curFileName));
end