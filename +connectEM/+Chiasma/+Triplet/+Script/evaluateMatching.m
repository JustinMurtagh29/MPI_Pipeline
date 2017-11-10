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

showPlots = true;

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

chiasmata = load(chiasmataFile);
axonFile = chiasmata.info.param.axonFile;
chiasmata = chiasmata.chiasmata;

axons = load(axonFile);
axons = axons.axons(axons.indBigAxons);

%% build chiasma table for triplets
chiasmaT = connectEM.Chiasma.Detect.buildTable(chiasmata);
chiasmaT(chiasmaT.nrExits ~= 3, :) = [];

%% calculate pair-wise distance
chiasmaT.posVx = cell2mat(arrayfun(@(a, c) ...
    chiasmata{a}.nodes(chiasmata{a}.ccCenterIdx(c), :), ...
    chiasmaT.aggloId, chiasmaT.chiasmaId, 'UniformOutput', false));
chiasmaT.posNm = chiasmaT.posVx .* param.raw.voxelSize;

% calculate pair-wise distance matrix
distMat = squareform(pdist(chiasmaT.posNm));

% set block diagonals to inifinity
blockRanges = reshape( ...
    1:size(chiasmaT, 1), [], 1);
blockRanges = accumarray( ...
    chiasmaT.aggloId, blockRanges, [], @(rows) {rows}, {[]});

for curRange = reshape(blockRanges, 1, [])
    distMat(curRange{1}, curRange{1}) = nan;
end

clear curRange blockRanges;

%% histogram over distance to closest candidate
if showPlots
minDistUm = min(distMat, [], 1) / 1E3;

figure;
histogram(minDistUm);
xlabel('Distance to closest triplet [µm]');
end

%% histogram over number of candidates within 2 µm
if showPlots
countWithinDist = sum(distMat < 2E3, 1);

figure;
histogram(countWithinDist);
xlabel('# triplets within 2 µm');
end

%% show examples with exactly one candidate within 2 µm
maskMat = (distMat < 2E3);
chiasmaIds = find(sum(maskMat, 1) == 1);
[partnerIds, ~] = find(maskMat(:, chiasmaIds));

rng(0);
randIds = randperm(numel(chiasmaIds));

if ~exist(outputDir, 'dir')
    mkdir(outputDir)
end

for curIdx = 1:20
    curIdOne = chiasmaIds(curIdx);
    curAggloIdOne = chiasmaT.aggloId(curIdOne);
    
    curIdTwo = partnerIds(curIdx);
    curAggloIdTwo = chiasmaT.aggloId(curIdTwo);
    
    curSkelOne = connectEM.Chiasma.Detect.buildSkeleton( ...
        axons(curAggloIdOne), chiasmata{curAggloIdOne}, chiasmaT.chiasmaId(curIdOne));
    curSkelTwo = connectEM.Chiasma.Detect.buildSkeleton( ...
        axons(curAggloIdTwo), chiasmata{curAggloIdTwo}, chiasmaT.chiasmaId(curIdTwo));
    curSkel = curSkelOne.mergeSkels(curSkelTwo);
    
    curSkel = Skeleton.setParams4Pipeline(curSkel, param);
    curFileName = sprintf('%02d_skel-pair.nml', curIdx);
    curSkel.write(fullfile(outputDir, curFileName));
end