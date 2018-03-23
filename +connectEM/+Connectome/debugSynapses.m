% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a_ax_spine_syn_clust.mat');
synFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v3_ax_spine_clustered.mat');
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_01_spine_attachment.mat');

outDir = '/home/amotta/Desktop';

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);
segPoints = Seg.Global.getSegToPointMap(param);

syn = load(synFile);
conn = load(connFile);

shAgglos = load(shFile, 'shAgglos');
shAgglos = shAgglos.shAgglos;

%% Absolute numbers
allSynCount = numel(syn.isSpineSyn);
allSpineSynCount = sum(syn.isSpineSyn);
allSpineSynFrac = mean(syn.isSpineSyn);

shLUT = Agglo.buildLUT(maxSegId, shAgglos);

% Reproduce spine flag
syn.mySpineIds = cellfun( ...
    @(ids) reshape(setdiff(shLUT(ids), 0), [], 1), ...
    syn.synapses.postsynId, 'UniformOutput', false);
syn.myIsSpine = ~cellfun(@isempty, syn.mySpineIds);
assert(isequal(syn.myIsSpine, syn.isSpineSyn));

allMultiSpineCount = sum( ...
    cellfun(@numel, syn.mySpineIds) > 1);

fprintf( ...
    '# synapses: %d\n', allSynCount);
fprintf( ...
    '# synapses onto spines: %d (%.0f %%)\n', ...
    allSpineSynCount, 100 * allSpineSynFrac);
fprintf( ...
    '# synapses onto multiple spines: %d\n', ...
    allMultiSpineCount);

%% Analyse spine heads
shCount = numel(shAgglos);

synShMat = [ ...
    repelem( ...
        transpose(1:allSynCount), ...
        cellfun(@numel, syn.mySpineIds)), ...
    cell2mat(syn.mySpineIds)];

[uniShIds, ~, uniSynIds] = unique(synShMat(:, 2));
uniSynIds = accumarray( ...
    uniSynIds, synShMat(:, 1), [], @(ids) {ids});
clear synShMat;

shT = table;
shT.id = reshape(1:shCount, [], 1);
shT.synIds(:) = {zeros(0, 1)};
shT.synIds(uniShIds) = uniSynIds;
clear uniShIds uniSynIds;

shWithSyn = ~cellfun(@isempty, shT.synIds);
shWithSynCount = sum(shWithSyn);
shWithSynFrac = mean(shWithSyn);

shWithMultiSynCount = sum( ...
    cellfun(@numel, shT.synIds) > 1);

fprintf( ...
    '# spine heads: %d\n', shCount);
fprintf( ...
    '# spine heads with synapse: %d (%.0f %%)\n', ...
    shWithSynCount, 100 * shWithSynFrac);
fprintf( ...
    '# spine heads with multiple synapses: %d\n', ...
    shWithMultiSynCount);

%% Look at spine heads with multiple synapses
shIds = shT.id(cellfun(@numel, shT.synIds) > 1);

% Select random subset
rng(0);
shIds = shIds(randperm(numel(shIds), 20));

% Generate skeleton representation
skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

for curIdx = 1:numel(shIds)
    curShId = shIds(curIdx);
    curShAgglo = shAgglos(curShId);

    curSynIds = shT.synIds{curShId};    
    curSynAgglos = syn.synapses(curSynIds, :);
    curSynAgglos = cellfun( ...
        @vertcat, ...
        curSynAgglos.presynId, ...
        curSynAgglos.postsynId, ...
        'UniformOutput', false);
    
    curNodes = cellfun( ...
        @(segIds) segPoints(segIds, :), ...
       [curShAgglo; curSynAgglos], ...
        'UniformOutput', false);
    skel = Skeleton.fromMST( ...
        curNodes, param.raw.voxelSize, skel);
    
    curPrefix = sprintf( ...
        '%0*d.', ceil(log10(1 + numel(shIds))), curIdx);
    curNames = [ ...
       {sprintf('Spine head #%d', curShId)}; ...
        arrayfun( ...
            @(id) sprintf('Synapse #%d', id), ...
            curSynIds, 'UniformOutput', false)];
    curNames = strcat(curPrefix, {' '}, curNames);
    skel.names((end - (numel(curNames) - 1)):end) = curNames;
end

skel.write(fullfile(outDir, 'multi-synaptic-spine-heads.nml'));
