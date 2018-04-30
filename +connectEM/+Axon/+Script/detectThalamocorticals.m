% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');

% These axons were manually identified based on the following tracings:
% * TC: https://webknossos.brain.mpg.de/annotations/Explorational/5ae5fff32700006943b3e2a8
% * CC: https://webknossos.brain.mpg.de/annotations/Explorational/5ae6c9c9270000c253b3f60b
tcAxonIds = [ ...
    15197, 4391, 185, 6758, 5761, ...
    4793, 9960, 5712, 7195, 2162];
ccAxonIds = [ ...
    20115, 9052, 19696, 21957, 5040, ...
    13742, 7553, 28847, 20707, 18837];

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

[connDir, connName] = fileparts(connFile);
interSynFile = fullfile(connDir, sprintf('%s_intersynapse.mat', connName));

conn = load(connFile);
interSyn = load(interSynFile);

syn = load(conn.info.param.synFile);
conn.axonMeta = connectEM.Axon.completeSynapseMeta(param, conn, syn);
    
%% calculate inter-synapse distances
% Previously we've calculate the synapse-to-synapse distances along the
% axon. To calculate the inter-synapse distances we now need to construct a
% tree-like representation of the axon. Let's do this by calculating the
% minimal spanning-tree over the synapses.
conn.axonMeta.pathLen = nan(size(conn.axonMeta.id));
conn.axonMeta.interSynDists = cell(size(conn.axonMeta.id));

for curIdx = 1:numel(interSyn.axonIds)
    curAxonId = interSyn.axonIds(curIdx);
    curPathLen = interSyn.axonPathLens(curIdx) / 1E3;
    
    curSynToSynDists = interSyn.synToSynDists{curIdx};
    curSynToSynDists = curSynToSynDists ./ 1E3;
    
    % NOTE(amotta): Zeros in the adjacency matrix are interpreted as
    % missing edge. This is a problem since synapse-to-synapse distance
    % zero occurs naturally when multiple synapses originate from the same
    % segment. Let's instead set these zeros to the smallest possible
    % non-zero value.
    curSynToSynDists(~curSynToSynDists) = eps;
    
    % claculate inter-synapse distances
    curInterSynDists = graphminspantree(sparse(curSynToSynDists));
    curInterSynDists = nonzeros(curInterSynDists);
    
    conn.axonMeta.pathLen(curAxonId) = curPathLen;
    conn.axonMeta.interSynDists{curAxonId} = curInterSynDists;
end

%% Plot synapse densities for GT axons
binEdges = linspace(0, 0.4, 21);

priSpineDens = ...
    conn.axonMeta.fullPriSpineSynCount ...
 ./ conn.axonMeta.pathLen;

fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [400, 400];

ax = axes(fig);
hold(ax, 'on');
axis(ax, 'square');

histogram( ...
    ax, priSpineDens(tcAxonIds), ...
    'BinEdges', binEdges, ...
    'DisplayStyle', 'stairs', ...
    'LineWidth', 2, ...
    'FaceAlpha', 1);
histogram( ...
    ax, priSpineDens(ccAxonIds), ...
    'BinEdges', binEdges, ...
    'DisplayStyle', 'stairs', ...
    'LineWidth', 2, ...
    'FaceAlpha', 1);

ax.TickDir = 'out';
xlim(ax, binEdges([1, end]));
yticks(ax, 0:ceil(ax.YLim(2)));

leg = legend(ax, ...
    'Thalamocortical', ...
    'Corticocortical', ...
    'Location', 'NorthEast');
leg.Box = 'off';

xlabel(ax, 'Spine synapse density (µm^{-1})');
ylabel(ax, 'Axons');

title( ...
    ax, {info.filename, info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
