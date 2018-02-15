% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a.mat');

minSynPre = 10;
spineFracThresh = 0.5;

% IDs of interneuron somata
% relative to all the 'Somata' agglomerates
inSomaIdsRel = [13, 56];
minInSomaSynCount = 1;
minSpineFracExc = 0.7;

info = Util.runInfo();

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

points = Seg.Global.getSegToPointMap(param);
sizes = Seg.Global.getSegToSizeMap(param);

conn = load(connFile);

%% find soma agglomerates
somaIds = find(conn.denMeta.targetClass == 'Somata');

%% soma synapse analysis
% show number of synapses per soma
fig = figure;
ax = axes(fig);

somaSynCount = conn.denMeta.synCount(somaIds);

histogram(ax, somaSynCount, 21);

xlabel('Synapses');
ylabel('Somata');

ax.XAxis.TickDirection = 'out';
ax.YAxis.TickDirection = 'out';

fig.Position(3:4) = [570, 350];

%% calculate per-soma position
somaAgglos = conn.dendrites(somaIds);
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

%% find axons which make synapses onto IN somata
inSomaConn = conn.connectome;

inSomaIds = find(conn.denMeta.targetClass == 'Somata');
inSomaIds = inSomaIds(inSomaIdsRel);

[~, inSomaConn.edges(:, 2)] = ...
    ismember(inSomaConn.edges(:, 2), inSomaIds);
inSomaConn(~inSomaConn.edges(:, 2), :)  = [];

inSomaAxonIds = find(conn.axonMeta.synCount >= minSynPre);
inSomaAxonIds = intersect(inSomaAxonIds, inSomaConn.edges(:, 1));

[~, inSomaConn.edges(:, 1)] = ...
    ismember(inSomaConn.edges(:, 1), inSomaAxonIds);
inSomaConn(~inSomaConn.edges(:, 1), :) = [];

% count synapses
inSomaConn.synCount = cellfun(@numel, inSomaConn.synIdx);

% count IN soma synapses per axon
axonInSomaSynCount = accumarray( ...
    inSomaConn.edges(:, 1), ...
    inSomaConn.synCount, ...
   [numel(inSomaAxonIds), 1]);

%% plot spine synapse fraction for these axons
spineSynFracAll = conn.axonMeta.spineSynCount ./ conn.axonMeta.synCount;
spineSynFrac = spineSynFracAll(inSomaAxonIds);

fig = figure();
ax = axes(fig);

histogram(ax, spineSynFrac, linspace(0, 1, 21));
xlim(ax, [0, 1]); xticks(ax, [0, 1]);
xlabel(ax, 'Spine synapse fraction');
ylabel(ax, {'Axons with'; '≥ 1 synapse onto IN somata'});

annotation( ...
    fig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
	'String', { ...
        'Spine synapse fraction of IN soma targeting axons'; ...
        info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');

fig.Position(3:4) = [678, 475];

%% plot specificities for these axons vs. all excitatory axons
excAxonMask = ...
    conn.axonMeta.synCount >= minSynPre ...
  & spineSynFracAll >= minSpineFracExc;
axonSpecsAll = conn.classConnectome ./ sum(conn.classConnectome, 2);

axonClasses = { ...
   {'Axons with'; '≥ 1 synapse onto IN somata'}, axonSpecsAll(inSomaAxonIds, :); ...
   {'Axons with'; sprintf('spine synapse fraction ≥ %.1f', minSpineFracExc)}, axonSpecsAll(excAxonMask, :)};
axonClassCount = size(axonClasses, 1);

fig = figure();

targetClassNames = conn.denClasses;
targetClassCount = numel(targetClassNames);

for curTargetClassIdx = 1:targetClassCount
    for curAxonClassIdx = 1:axonClassCount
        curAx = subplot( ...
            axonClassCount, targetClassCount, ...
            curTargetClassIdx + (curAxonClassIdx - 1) * targetClassCount);
        
        curAxonSpecs = axonClasses{curAxonClassIdx, 2};
        curAxonSpecs = curAxonSpecs(:, curTargetClassIdx);
        histogram(curAx, curAxonSpecs, linspace(0, 1, 21));
        
        curAx.TickDir = 'out';
        
        xlim(curAx, [0, 1]);
        xticks(curAx, [0, 1]);
        xticklabels(curAx, {});
        
        if curAxonClassIdx == axonClassCount
            xticklabels(curAx, {'0', '1'});
            xlabel(curAx, sprintf('S(%s)', ...
                targetClassNames{curTargetClassIdx}));
        end
        
        if curTargetClassIdx == 1
            curAxonClassName = axonClasses{curAxonClassIdx, 1};
            ylabel(curAx, curAxonClassName);
        end
    end
end

annotation( ...
    fig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
	'String', { ...
        'IN soma targeting axons versus all spine-targeting axons'; ...
        info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');

fig.Position(3:4) = [1758, 631];