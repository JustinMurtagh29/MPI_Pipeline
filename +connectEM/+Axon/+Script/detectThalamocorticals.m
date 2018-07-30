% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', ...
    'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

% Ground truth annotations
oldConnFile = fullfile(rootDir, 'connectomeState', ...
    'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');

% These axons were manually identified based on the following tracings:
% * TC: https://webknossos.brain.mpg.de/annotations/Explorational/5ae5fff32700006943b3e2a8
% * CC: https://webknossos.brain.mpg.de/annotations/Explorational/5ae6c9c9270000c253b3f60b
oldTcAxonIds = [ ...
    15197, 4391, 185, 6758, 5761, ...
    4793, 9960, 5712, 7195, 2162];
oldCcAxonIds = [ ...
    20115, 9052, 19696, 21957, 5040, ...
    13742, 7553, 28847, 20707, 18837];

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

[conn, syn, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);

[connDir, connName] = fileparts(connFile);
interSynFile = sprintf('%s_intersynapse_v2.mat', connName);
interSynFile = fullfile(connDir, interSynFile);
interSyn = load(interSynFile);

%% Translate ground truth annotations to latest connectome
clear cur*;
curConn = load(oldConnFile);
curMaxSegId = Seg.Global.getMaxSegId(param);

curNewLUT = Agglo.buildLUT(curMaxSegId, conn.axons);
curFindInNew = @(segIds) mode(nonzeros(curNewLUT(segIds)));

tcAxonIds = cellfun(curFindInNew, curConn.axons(oldTcAxonIds));
tcAxonIds = reshape(tcAxonIds, size(oldTcAxonIds));

ccAxonIds = cellfun(curFindInNew, curConn.axons(oldCcAxonIds));
ccAxonIds = reshape(ccAxonIds, size(oldCcAxonIds));
clear cur*;

%% Calculate features for TC detection
synIds = connectEM.Axon.getSynapses(conn, syn);

synIsSpine = cellfun( ...
    @(ids) syn.synapses.type(ids) == 'PrimarySpine', ...
    synIds, 'UniformOutput', false);

boutonIds = ...
    connectEM.Axon.clusterSynapsesIntoBoutons( ...
        synIds, interSyn, 'cutoffDist', 2433);

conn.axonMeta.fullPriSpinesPerBouton = cellfun( ...
    @(boutonIds, isSpine) mean(accumarray(boutonIds, isSpine)), ...
    boutonIds, synIsSpine);

conn.axonMeta.fullPriSpineSynDens = ...
    conn.axonMeta.fullPriSpineSynCount ...
 ./ conn.axonMeta.pathLen;

%% Plot synapse densities for GT axons
thresh = 0.19;
excAxonIds = axonClasses(1).axonIds;
binEdges = linspace(0, 0.4, 21);

fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [425, 400];

ax = axes(fig);
colors = ax.ColorOrder;

hold(ax, 'on');
axis(ax, 'square');

yyaxis(ax, 'right');
tcHist = histogram(ax, conn.axonMeta.fullPriSpineSynDens(tcAxonIds));
ccHist = histogram(ax, conn.axonMeta.fullPriSpineSynDens(ccAxonIds));
ylim(ax, [0, 1]);
ylabel(ax, 'Probability');

yyaxis(ax, 'left');
excHist = histogram(ax, conn.axonMeta.fullPriSpineSynDens(excAxonIds));

allHists = [tcHist, ccHist, excHist];
[allHists.DisplayStyle] = deal('stairs');
[allHists.BinEdges] = deal(binEdges);
[allHists.LineWidth] = deal(2);
[allHists.FaceAlpha] = deal(1);

[allHists(1:2).Normalization] = deal('probability');
allHists(1).EdgeColor = colors(1, :);
allHists(2).EdgeColor = colors(2, :);
allHists(3).EdgeColor = zeros(1, 3);

plot( ...
    ax, [thresh, thresh], ax.YLim, ...
    'Color', 'black', 'LineStyle', '--');

ax.TickDir = 'out';
ax.YColor = [0, 0, 0];
xlim(ax, binEdges([1, end]));
xlabel(ax, 'Spine synapse density (µm^{-1})');
ylabel(ax, 'Excitatory axons');

leg = legend( ...
    [tcHist, ccHist], ...
    'Thalamocortical', ...
    'Corticocortical', ...
    'Location', 'NorthEast');
leg.Box = 'off';

title( ...
    ax, {info.filename, info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
