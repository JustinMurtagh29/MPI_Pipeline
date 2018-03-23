% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
param = struct;
param.saveFolder = '/gaba/u/mberning/results/pipeline/20170217_ROI';

connFile = fullfile(param.saveFolder, 'connectomeState', 'connectome_axons_18_a_ax_spine_syn_clust.mat');
synFile = fullfile(param.saveFolder, 'connectomeState', 'SynapseAgglos_v3_ax_spine_clustered.mat');
shFile = fullfile(param.saveFolder, 'aggloState', 'dendrites_wholeCells_01_spine_attachment.mat');

info = Util.runInfo();

%% Loading data
maxSegId = Seg.Global.getMaxSegId(param);

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
