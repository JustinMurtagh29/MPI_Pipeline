% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
axonFile = fullfile(rootDir, 'aggloState', 'axons_09_a.mat');

chiasmaDir = fullfile( ...
    rootDir, 'chiasmataSplitting', ...
    '20171102T150211-ortho-mode-on-axons-9');
chiasmaFile = fullfile(chiasmaDir, '20171102T150211_chiasmata.mat');

% chiasmata
chiParam = struct;
chiParam.sphereRadiusOuter = 10000; % in nm
chiParam.sphereRadiusInner = 1000; % in nm
chiParam.minNodeDist = 2000; % in nm
chiParam.clusterSize = 2000; % in nm

curDateStr = datestr(now, 30);
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

%% build NMLs
rng(0);
randIds = randperm(size(skelData, 1));
taskDefs = struct([]);

nmlDir = sprintf('%s_queries', curDateStr);
nmlDir = fullfile(chiasmaDir, nmlDir);
if ~exist(nmlDir, 'dir'); mkdir(nmlDir); end

for curIdx = 1:numel(randIds)
    curRow = randIds(curIdx);
    
    curSkelData = skelData(curRow, :);
   [curSkel, curTaskDef] = connectEM.Chiasma.Ortho.buildQuery( ...
        param, axons(curSkelData.axonId), curSkelData.nodeId);
    
    curNmlFile = sprintf( ...
        '%d_axon-%d_node-%d.nml', curIdx, ...
        curSkelData.axonId, curSkelData.nodeId);
    curSkel.write(fullfile(nmlDir, curNmlFile));
    
    curTaskDef.axonId = curSkelData.axonId;
    curTaskDef.nmlFile = curNmlFile;
    
    if curIdx == 1; taskDefs = curTaskDef; end
    taskDefs(curIdx) = curTaskDef;
end

taskDefs = orderfields(taskDefs);
taskDefs = reshape(taskDefs, [], 1);

%% save query data
out = struct;
out.info = info;
out.taskDefs = taskDefs;

outFile = sprintf('%s_tasks.mat', curDateStr);
outFile = fullfile(chiasmaDir, outFile);

Util.saveStruct(outFile, out);
system(sprintf('chmod a-w "%s"', outFile));