% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

connFile = fullfile( ...
    rootDir, 'connectomeState', ...
    'connectome_axons_18_a_with_den_meta.mat');

minSynPre = 10;
minSynPost = 0;
spineFracThresh = 0.5;

%% loading data
conn = load(connFile);

%% soma synapse analysis
% show number of synapses per soma

fig = figure;
ax = axes(fig);

somaMask = (conn.denMeta.targetClass == 'Somata');
somaSynCount = conn.denMeta.synCount(somaMask);

histogram(ax, somaSynCount, 21);

xlabel('Synapses');
ylabel('Somata');

ax.XAxis.TickDirection = 'out';
ax.YAxis.TickDirection = 'out';

fig.Position(3:4) = [570, 350];