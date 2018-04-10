% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_01_spine_attachment.mat');
synFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v3_ax_spine_clustered.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a_ax_spine_syn_clust.mat');

% As in Benedikt's +L4/updateParamsToNewestFiles.m
% Commit hash `590d8538d65463151d43491e2446e25ca11dd5f6`
graphFile = fullfile(rootDir, 'graphNew.mat');
synScoreFile = fullfile(rootDir, 'globalSynScores.mat');

minSynCount = 10;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

% Prepare for load synapse scores
param.svg.graphFile = graphFile;
param.svg.synScoreFile = synScoreFile;
graph = Seg.IO.loadGraph(param, false);

maxSegId = Seg.Global.getMaxSegId(param);
segPoints = Seg.Global.getSegToPointMap(param);

shAgglos = load(shFile, 'shAgglos');
shAgglos = shAgglos.shAgglos;

syn = load(synFile);
conn = load(connFile);

%% Classify synapses
somaAgglos = conn.denMeta.targetClass == 'Somata';
somaAgglos = conn.dendrites(somaAgglos);

syn.synapses.type = ...
    connectEM.Synapse.classifyType( ...
        param, syn.synapses, graph.synScores, ...
        shAgglos, somaAgglos, conn.axons);

%% Collect synapses
shLUT = (Agglo.buildLUT(maxSegId, shAgglos) ~= 0);

syn.synapses.id = reshape( ...
    1:size(syn.synapses, 1), [], 1);
syn.synapses.ontoSpine = cellfun( ...
    @(segIds) any(shLUT(segIds)), ...
    syn.synapses.postsynId);
clear shLUT;

% This time we count a synapse even when none of the postsynaptic segments
% is in a dendrite agglomerates. This should give a more accurate picture
% of the axonal outputs.
axonLUT = Agglo.buildLUT(maxSegId, conn.axons);

synapses = syn.synapses;
synapses.axonId = cellfun( ...
    @(segIds) setdiff(axonLUT(segIds), 0), ...
    synapses.presynId, 'UniformOutput', false);
clear axonLUT;

% Remove synapses which have
% * no presynaptic axon at all
% * multiple axons on presynaptic side (?!)
synapses(~cellfun(@isscalar, synapses.axonId), :) = [];
synapses.axonId = cell2mat(synapses.axonId);

axonMeta = conn.axonMeta;
axonMeta.fullSynCount = accumarray( ...
    synapses.axonId, 1, size(axonMeta.id));
axonMeta.fullSpineSynCount = accumarray( ...
    synapses.axonId, synapses.ontoSpine, size(axonMeta.id));
axonMeta.fullPriSpineSynCount = accumarray( ...
    synapses.axonId, synapses.type == 'PrimarySpine', size(axonMeta.id));

axonMeta.synIds = accumarray( ...
    conn.connectome.edges(:, 1), ...
    transpose(1:size(conn.connectome, 1)), size(axonMeta.id), ...
    @(r) {cell2mat(conn.connectome.synIdx(r))}, {zeros(0, 1)});
axonMeta.fullSynIds = accumarray( ...
    synapses.axonId, synapses.id, size(axonMeta.id), ...
    @(synIds) {synIds}, {zeros(0, 1)});

% Sanity checks
assert(all(cellfun(@numel, axonMeta.synIds) == axonMeta.synCount));
assert(all(cellfun(@numel, axonMeta.fullSynIds) == axonMeta.fullSynCount));

axonMeta.spineSynFrac = ...
    axonMeta.spineSynCount ...
 ./ axonMeta.synCount;
axonMeta.fullSpineSynFrac = ...
    axonMeta.fullSpineSynCount ...
 ./ axonMeta.fullSynCount;
axonMeta.fullPriSpineSynFrac = ...
    axonMeta.fullPriSpineSynCount ...
 ./ axonMeta.fullSynCount;

clear synapses;

% Remove axons with too few synapses
axonMeta(axonMeta.synCount < minSynCount, :) = [];

%% Plot spine synapse fractions
binEdges = linspace(0, 1, 21);

fig = figure();
fig.Color = 'white';

ax = axes(fig);
axis(ax, 'square')
hold(ax, 'on');
histogram(ax, ...
    axonMeta.fullPriSpineSynFrac, binEdges, ...
    'DisplayStyle', 'stairs', 'LineWidth', 2);

ax.TickDir = 'out';
xlim(ax, binEdges([1, end]));
xlabel(ax, 'Spine synapse fraction');
ylabel(ax, 'Axons');

title(ax, ...
   {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
