% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
chiasmataFile = fullfile( ...
    rootDir, 'chiasmataSplitting', ...
    '20171104T181213-on-axons-10a-plus-10kE3a', ...
    '20171104T184018_chiasmata.mat');

outputDir = '/home/amotta/Desktop/random-chiasmata';

%% load data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

data = load(chiasmataFile);
axonFile = data.info.param.axonFile;
chiasmata = data.chiasmata;
clear data;

axons = load(axonFile);
axonIds = find(axons.indBigAxons);
axons = axons.axons(axonIds);

%% select chiasmata
chiasma = table;
chiasma.axonId = repelem((1:numel(chiasmata))', ...
    cellfun(@(c) numel(c.ccCenterIdx), chiasmata));
chiasma.chiasmaId = cell2mat(cellfun(@(c) ...
    (1:numel(c.ccCenterIdx))', chiasmata, 'UniformOutput', false));

% restrict to 4-fold chiasmata
chiasma.nrExits = arrayfun(@(a, c) ...
    chiasmata{a}.nrExits(chiasmata{a}.ccCenterIdx(c)), ...
    chiasma.axonId, chiasma.chiasmaId);
chiasma(chiasma.nrExits ~= 4, :) = [];
chiasma.nrExits = [];

% remove solved chiasmata
chiasma.isSolved = arrayfun(@(a, c) ...
    axons(a).solvedChiasma(chiasmata{a}.ccCenterIdx(c)), ...
    chiasma.axonId, chiasma.chiasmaId);
chiasma(chiasma.isSolved, :) = [];
chiasma.isSolved = [];

%% export axons
rng(0);
randIds = randperm(size(chiasma, 1));

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

for curIdx = 1:20
    curChiasma = chiasma(randIds(curIdx), :);
    
    curSkel = connectEM.Chiasma.Detect.buildSkeleton( ...
        axons(curChiasma.axonId), chiasmata{curChiasma.axonId});
    curSkel = Skeleton.setParams4Pipeline(curSkel, param);
    
    curSkelFile = sprintf( ...
        '%03d_axon-%d_chiasma-%d.nml', ...
        curIdx, curChiasma.axonId, curChiasma.chiasmaId);
    curSkel.write(fullfile(outputDir, curSkelFile));
end