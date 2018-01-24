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
info = Util.runInfo();

%% loading data
dendData = load(dendFile);
sdData = load(sdFile);
adData = load(adFile);
aisData = load(aisFile);

%% sanity checks
idxBig = find(dendData.indBigDends(:));
numBigDends = sum(dendData.indBigDends);

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
