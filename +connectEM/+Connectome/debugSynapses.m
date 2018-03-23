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

%% Build spine head look-up table
shLUT = Agglo.buildLUT(maxSegId, shAgglos);

%% Absolute numbers
allSynCount = numel(syn.isSpineSyn);
allSpineSynCount = sum(syn.isSpineSyn);
allSpineSynFrac = mean(syn.isSpineSyn);

% Reproduce spine flag
syn.myIsSpine = cellfun( ...
    @(ids) any(shLUT(ids)), syn.synapses.postsynId);
assert(isequal(syn.myIsSpine, syn.isSpineSyn));