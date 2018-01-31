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
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

points = Seg.Global.getSegToPointMap(param);
sizes = Seg.Global.getSegToSizeMap(param);

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

%% calculate per-soma position
somaAgglos = conn.dendrites(somaMask);
somaSize = cellfun(@(segIds) ...
    sum(sizes(segIds)), somaAgglos);
somaPos = cellfun(@(segIds) ...
    sum(sizes(segIds) .* points(segIds, :), 1), ...
    somaAgglos, 'UniformOutput', false);
somaPos = ceil(cell2mat(somaPos) ./ somaSize);

%% export somata to webKNOSSOS
[~, somaIds] = sort(somaSynCount, 'descend');

numDigits = ceil(log10(1 + numel(somaSynCount)));
treeNames = arrayfun(@(r, i, n) sprintf( ...
    '%0*d. Soma %d (%d synapses)', numDigits, r, i, n), ...
    reshape(1:numel(somaSynCount), [], 1), ...
    somaIds, somaSynCount(somaIds, :), ...
    'UniformOutput', false);

skel = skeleton();
skel = skel.addNodesAsTrees( ...
    somaPos(somaIds, :), treeNames);
skel = Skeleton.setParams4Pipeline(skel, param);
skel.write('/home/amotta/Desktop/somata-with-synapse-count.nml');