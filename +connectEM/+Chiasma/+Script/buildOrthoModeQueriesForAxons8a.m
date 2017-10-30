% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
axonFile = fullfile(rootDir, 'aggloState', 'axons_08_b.mat');
chiasmaFile = '/mnt/mpibr/data/Personal/mottaa/L4/2017-10-30-Chiasmata-Detection-On-Axons-8a/chiasmata.mat';
outputDir = '/home/amotta/Desktop/ortho-chiasma-queries';

% chiasmata
chiParam = struct;
chiParam.sphereRadiusOuter = 10000; % in nm
chiParam.sphereRadiusInner = 1000; % in nm
chiParam.minNodeDist = 2000; % in nm
chiParam.clusterSize = 2000; % in nm

info = Util.runInfo();

%% load data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

axons = load(axonFile, 'axons', 'indBigAxons');
axonIds = find(axons.indBigAxons);
axons = axons.axons(axonIds);

chiasmata = load(chiasmaFile, 'chiasmata');
chiasmata = chiasmata.chiasmata;

% sanity checks
assert(numel(axons) == numel(chiasmata));

%% add chiasma parameters
chiParam.voxelSize = param.raw.voxelSize;

chiParam = cat(2, fieldnames(chiParam), struct2cell(chiParam));
chiParam = reshape(transpose(chiParam), 1, []);

param = Util.modifyStruct(param, chiParam{:});

%% 
nodeIds = cellfun( ...
    @(s) s.ccCenterIdx(:), ...
    chiasmata, 'UniformOutput', false);

skelData = table;
skelData.axonId = repelem( ...
    (1:numel(nodeIds))', cellfun(@numel, nodeIds));
skelData.nodeId = cell2mat(nodeIds(:));
clear nodeIds;

% skip solved chiasmata
skelData.solved = arrayfun( ...
    @(aIdx, nIdx) axons(aIdx).solvedChiasma(nIdx), ...
    skelData.axonId, skelData.nodeId);
skelData(skelData.solved, :) = [];

%%
rng(0);
randIds = randperm(size(skelData, 1));

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

for curIdx = 1:numel(randIds)
    curRow = randIds(curIdx);
    
    curSkelData = skelData(curRow, :);
    curSkel = connectEM.Chiasma.buildOrthoModeQuery( ...
        param, axons(curSkelData.axonId), curSkelData.nodeId);
    
    curSkelFile = sprintf( ...
        '%d_axon-%d_node-%d.nml', curIdx, ...
        curSkelData.axonId, curSkelData.nodeId);
    curSkel.write(fullfile(outputDir, curSkelFile));
end