% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a_with_den_meta.mat');

%% loading data
% NOTE(amotta): Synapses sizes are contained in the `contactArea` field of 
% `conn.connectomeMeta`. Each cell contains the synapses sizes of the
% correponding entries in `conn.connectome`.
conn = load(connFile);

%% find doubly coupled neurites
conn.connectome.synCount = ...
    cellfun(@numel, conn.connectome.synIdx);

prePostPairMask = (conn.connectome.synCount == 2);
nrPrePostPairs = sum(prePostPairMask);

jointPreIds = unique(conn.connectome.edges(prePostPairMask, 1));
jointPostIds = unique(conn.connectome.edges(prePostPairMask, 2));

%% plot correlation for pre- & post-coupled synapses
prePostPairs = conn.connectome.edges(prePostPairMask, :);
prePostPairSynAreas = conn.connectomeMeta.contactArea(prePostPairMask);

% sort areas within
prePostPairSynAreas = cellfun( ...
    @(areas) reshape(sort(areas), 1, []), ...
    prePostPairSynAreas, 'UniformOutput', false);

fig = figure();
ax = axes(fig);

data = cell2mat(prePostPairSynAreas);
dataMax = max(data(:));

scatter(ax, data(:, 2), data(:, 1), 18, '+');
xlim(ax, [0, dataMax]);
ylim(ax, [0, dataMax]);