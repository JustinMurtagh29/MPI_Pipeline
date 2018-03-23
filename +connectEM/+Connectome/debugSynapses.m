% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
param = struct;
param.saveFolder = '/gaba/u/mberning/results/pipeline/20170217_ROI';

connName = 'connectome_axons_18_a_ax_spine_syn_clust';
synFile = fullfile(param.saveFolder, 'connectomeState', 'SynapseAgglos_v3_ax_spine_clustered.mat');
shFile = fullfile(param.saveFolder, 'aggloState', 'dendrites_wholeCells_01_spine_attachment.mat');

info = Util.runInfo();

%% Loading data
maxSegId = Seg.Global.getMaxSegId(param);

syn = load(synFile);
conn = connectEM.Connectome.load(param, connName);

shAgglos = load(shFile, 'shAgglos');
shAgglos = shAgglos.shAgglos;

%% Absolute numbers
allSynCount = numel(syn.isSpineSyn);
allSpineSynCount = sum(syn.isSpineSyn);
allSpineSynFrac = mean(syn.isSpineSyn);

shLUT = Agglo.buildLUT(maxSegId, shAgglos);

% Reproduce spine flag
syn.myIsSpine = cellfun( ...
    @(ids) any(shLUT(ids)), syn.synapses.postsynId);
assert(isequal(syn.myIsSpine, syn.isSpineSyn));

fprintf( ...
    '# synapses: %d\n', allSynCount);
fprintf( ...
    '# synapses onto spines: %d (%.0f %%)\n', ...
    allSpineSynCount, 100 * allSpineSynFrac);

%% Analyse spine heads
shCount = numel(shAgglos);

postSynSegIds = cell2mat(syn.synapses.postsynId);
postSynSegLUT = logical(Agglo.buildLUT(maxSegId, {postSynSegIds}));

shWithSyn = cellfun( ...
    @(ids) any(postSynSegLUT(ids)), shAgglos);
shWithSynCount = sum(shWithSyn);
shWithSynFrac = mean(shWithSyn);

fprintf( ...
    '# spine heads: %d\n', shCount);
fprintf( ...
    '# spine heads with synapse: %d (%.0f %%)\n', ...
    shWithSynCount, 100 * shWithSynFrac);