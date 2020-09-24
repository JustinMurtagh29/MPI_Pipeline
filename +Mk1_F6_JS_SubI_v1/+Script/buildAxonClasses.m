% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet/';
connFile = fullfile(rootDir, 'connectome', 'Connectome_20191227T220548-results_20191227T220548-results-auto-spines-v3_SynapseAgglomerates--20191227T220548-results--20191227T220548-results-auto-spines-v3--v1.mat');

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

conn = load(connFile);

synFile = conn.info.param.synFile;
syn = load(synFile);

curTypeFile = strrep(synFile, '.mat', '__types_v1.mat');
syn.synapses.type = Util.load(curTypeFile, 'types');

%% Build axon types
clear cur*;

synT = connectEM.Connectome.buildSynapseTable(conn, syn);
synT.type = syn.synapses.type(synT.id);

axonMeta = conn.axonMeta;

% NOTE(amotta): Remove "old" synapse counts
curVars = axonMeta.Properties.VariableNames;
curVars = setdiff(curVars, {'synCount', 'spineSynCount'}, 'stable');
axonMeta = axonMeta(:, curVars);

axonMeta.synCount = accumarray( ...
    synT.preAggloId, 1, [height(axonMeta), 1]);
axonMeta.priSpineSynCount = accumarray( ...
    synT.preAggloId, synT.type == 'PrimarySpine', [height(axonMeta), 1]);

axonMeta.priSpineSynFrac = ...
    axonMeta.priSpineSynCount ./ axonMeta.synCount;

%% Plot results
clear cur*;
curBinEdges = linspace(0, 1, 11);
curLeg = {};

fig = figure();
ax = axes(fig);
hold(ax, 'on');

curMask = axonMeta.synCount >= 5;
curLeg{end + 1} = '≥ 5 synapses';
histogram(axonMeta.priSpineSynFrac(curMask), curBinEdges);

curMask = axonMeta.synCount >= 10;
curLeg{end + 1} = '≥ 10 synapses';
histogram(axonMeta.priSpineSynFrac(curMask), curBinEdges);

curMask = axonMeta.synCount >= 10 & axonMeta.synCount <= 25;
curLeg{end + 1} = '10 - 25 synapses';
histogram(axonMeta.priSpineSynFrac(curMask), curBinEdges);

curHists = flip(ax.Children);
set(curHists, 'Normalization', 'probability');

% NOTE(amotta): Add axon counts to legend
curCounts = cellfun(@numel, get(curHists, {'Data'}));
curCounts = reshape(num2cell(curCounts), size(curLeg));

curLeg = cellfun( ...
    @(s, n) sprintf('%s (n = %d)', s, n), ...
    curLeg, curCounts, 'UniformOutput', false);

curLeg = legend(curHists, curLeg, 'Location', 'NorthWest');

xlabel(ax, 'Fraction of primary spine output synapses');
ylabel(ax, 'Fraction of axons');

fig.Position(3:4) = [440, 310];
connectEM.Figure.config(fig, info);

saveas(gcf, fullfile(rootDir, 'connectome','figures','buildAxonClasses.png'))
close all
