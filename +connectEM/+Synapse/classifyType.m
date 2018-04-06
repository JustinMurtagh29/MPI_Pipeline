% Classifies synapses into
% * soma synapses
% * shaft synapses
% * primary spine synapse
% * secondary spine synapse
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_01_spine_attachment.mat');
somaFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a_ax_spine_syn_clust.mat');
synFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v3_ax_spine_clustered.mat');

% As in Benedikt's +L4/updateParamsToNewestFiles.m
% Commit hash `590d8538d65463151d43491e2446e25ca11dd5f6`
graphFile = fullfile(rootDir, 'graphNew.mat');
synScoreFile = fullfile(rootDir, 'globalSynScores.mat');

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

% Prepare for load synapse scores
param.svg.graphFile = graphFile;
param.svg.synScoreFile = synScoreFile;

% Synapse scores

maxSegId = Seg.Global.getMaxSegId(param);
graph = Seg.IO.loadGraph(param, false);

%%

% Spine heads
shAgglos = load(shFile, 'shAgglos');
shAgglos = shAgglos.shAgglos;

% Somata
somaAgglos = load(somaFile, 'dendrites', 'denMeta');
somaAgglos = somaAgglos.dendrites( ...
    somaAgglos.denMeta.targetClass == 'Somata');

% Synapses
synapses = load(synFile, 'synapses');
synapses = synapses.synapses;

%%
% Building look-up tables
% Soma class dominates over spine heads
somaLUT = Agglo.buildLUT(maxSegId, somaAgglos);
shLUT = Agglo.buildLUT(maxSegId, shAgglos);
shLUT(somaLUT ~= 0) = 0;

synapses.id = reshape( ...
    1:size(synapses, 1), [], 1);
synapses.synScores = cellfun( ...
    @(ids) max(graph.synScores(ids, :), [], 2), ...
    synapses.edgeIdx, 'UniformOutput', false);

% Sanity check
% Synapse threshold was -1.67
assert(all(cellfun(@min, synapses.synScores) > -1.67));
synapses.maxSynScore = cellfun(@max, synapses.synScores);

synapses.somaId = cellfun( ...
    @(segIds) max(somaLUT(segIds)), ...
    synapses.postsynId);
synapses.shId = cellfun( ...
    @(segIds) max(shLUT(segIds)), ...
    synapses.postsynId);

%% Separate primary from secondary spine synapses
shSynIds = accumarray( ...
    1 + synapses.shId, synapses.id, ...
   [1 + numel(shAgglos), 1], @(ids) {ids}, {zeros(0, 1)});
shSynIds = shSynIds(2:end);

% Sort synapse by SynEM scores
[~, sortedIds] = cellfun( ...
    @(synIds) sort(synapses.maxSynScore(synIds)), ...
    shSynIds, 'UniformOutput', false);
shSynIds = cellfun( ...
    @(synIds, sortIds) synIds(sortIds), ...
    shSynIds, sortedIds, 'UniformOutput', false);
clear sortedIds;

shT = table;
shT.id = reshape( ...
    1:numel(shAgglos), [], 1);

mask = ~cellfun(@isempty, shSynIds);

shT.priSynId(:) = 0;
shT.priSynId(mask) = cellfun( ...
    @(synIds) synIds(end), shSynIds(mask));

shT.secSynIds(:) = {zeros(0, 1)};
shT.secSynIds(mask) = cellfun( ...
    @(synIds) synIds(1:(end - 1)), ...
    shSynIds(mask), 'UniformOutput', false);

clear mask;

%% Look at examples in webKNOSSOS