% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
dendFile = fullfile( ...
    rootDir, 'aggloState', 'dendrites_wholeCells_01_spine_attachment.mat');

sdFile  = '/tmpscratch/sahilloo/L4/dataPostSyn/dendritesSmoothState.mat';
adFile  = '/tmpscratch/sahilloo/L4/dataPostSyn/dendritesADState.mat';
aisFile = '/tmpscratch/sahilloo/L4/dataPostSyn/dendritesAISState.mat';

% set to export NML files
nmlDir = '';

info = Util.runInfo();

%% loading data
dendData = load(dendFile);
sdData = load(sdFile);
adData = load(adFile);
aisData = load(aisFile);

%% sanity checks
idxBig = find(dendData.indBigDends(:));
numBigDends = sum(dendData.indBigDends);
numWholeCells = numel(dendData.indWholeCells);

fprintf('Overall:\n');
fprintf('  # large dendrites: %d\n', numBigDends);
fprintf('  # whole cells: %d\n', numWholeCells);
fprintf('\n');

% for smooth dendrites
assert(numBigDends == numel(sdData.idxSmooth));
fprintf('Smooth dendrites\n');
fprintf('  # found: %d\n', sum(sdData.idxSmooth));
fprintf('\n');

% for apical dendrites
assert(all(adData.idxAD <= numBigDends));
assert(all(adData.idxAD > 0));
fprintf('Apical dendrites\n');
fprintf('  # found: %d\n', numel(adData.idxAD))
fprintf('  # also in smooth dendrites: %d\n', ...
    sum(sdData.idxSmooth(adData.idxAD)));
fprintf('\n');

% for axon initial segments
assert(all(aisData.idxAIS <= numBigDends));
assert(all(aisData.idxAIS > 0));
fprintf('Axon initial segments\n');
fprintf('  # found: %d\n', numel(aisData.idxAIS));
fprintf('  # also in smooth dendrites: %d\n', ...
    sum(sdData.idxSmooth(aisData.idxAIS)));
fprintf('  # also in apical dendrites: %d\n', ...
    numel(intersect(adData.idxAD, aisData.idxAIS)));
fprintf('\n');

%%
targetClass = zeros(size(dendData.dendAgglos));
targetClass(idxBig)                   = 1;
targetClass(idxBig(sdData.idxSmooth)) = 2;
targetClass(idxBig(adData.idxAD))     = 3;
targetClass(idxBig(aisData.idxAIS))   = 4;
targetClass(dendData.indWholeCells)   = 5;

targetLabels = { ...
    'Ignore', 'OtherDendrite', 'SmoothDendrite', ...
    'ApicalDendrite', 'AxonInitialSegment', 'WholeCell'};
targetClass = categorical( ...
    targetClass, 0:5, targetLabels);

%% export examples
if ~isempty(nmlDir)
    mkdir(nmlDir);
    
    param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
    param = param.p;
    
    voxelSize = param.raw.voxelSize;
    points = Seg.Global.getSegToPointMap(param);
    
    for curLabelIdx = 3:numel(targetLabels)
        curLabel = targetLabels{curLabelIdx};
        curAggloIds = find(targetClass(:) == curLabel);
        curAgglos = dendData.dendAgglos(curAggloIds);
        
        curNodes = cellfun( ...
            @(segIds) points(segIds, :), ...
            curAgglos, 'UniformOutput', false);
        
        curSkel = Skeleton.fromMST(curNodes, voxelSize);
        curSkel.names = arrayfun( ...
            @(i) sprintf('Dendrite #%d', i), ...
            curAggloIds, 'UniformOutput', false);
        
        curNmlFile = strcat(lower(curLabel), '.nml');
        curNmlFile = fullfile(nmlDir, curNmlFile);
        
        curSkel = Skeleton.setParams4Pipeline(curSkel, param);
        curSkel.write(curNmlFile);
    end
end
