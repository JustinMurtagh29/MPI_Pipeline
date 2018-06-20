% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');

targetClasses = { ...
    'Somata', 'SO'; ...
    'ProximalDendrite', 'PD'; ...
    'ApicalDendrite', 'AD'; ...
    'SmoothDendrite', 'SD'; ...
    'AxonInitialSegment', 'AIS'};
targetLabels = targetClasses(:, 2);
targetClasses = targetClasses(:, 1);

minSynPre = 10;
info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, ~, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);

%% Preparing data
axonClasses(1).tag = 'Exc';
axonClasses(2).tag = 'Inh';
axonClasses(3).tag = 'TC';
axonClasses(4).tag = 'CC';

[conn, axonClasses] = ...
    connectEM.Connectome.prepareForSpecificityAnalysis( ...
        conn, axonClasses, 'minSynPre', minSynPre);
axonClasses = ...
    connectEM.Connectome.buildAxonSpecificityClasses(conn, axonClasses);

axonSynCounts = nan(size(conn.axons));
axonSynCounts(conn.axonMeta.id) = conn.axonMeta.synCount;

%% Plotting
plotAxonClassIds = [4, 3, 2];
numTargetClasses = numel(targetClasses);

fig = figure();
fig.Color = 'white';

for curAxonClassIdx = 1:numel(plotAxonClassIds)
    curAxonClassId = plotAxonClassIds(curAxonClassIdx);
    curAxonClass = axonClasses(curAxonClassId);
    
    for curTargetClassId = 1:numTargetClasses
        curTargetClass = targetClasses{curTargetClassId};
        curSpecAxonIds = curAxonClass.specs.(curTargetClass).axonIds;
        curSynCounts = axonSynCounts(curSpecAxonIds);
        curBinEdges = 0:max(curSynCounts);
        
        curAx = subplot( ...
            numel(plotAxonClassIds), numTargetClasses, ...
            (curAxonClassIdx - 1) * numTargetClasses + curTargetClassId);
        histogram(curAx, ...
            curSynCounts, ...
            'BinEdges', curBinEdges, ...
            'FaceColor', curAx.ColorOrder(1, :), ...
            'EdgeColor', 'none', ...
            'FaceAlpha', 1);
        
        curAx.Box = 'off';
        curAx.TickDir = 'out';
    end
    
    curAxes = fig.Children(1:numTargetClasses);
    curMaxX = max(arrayfun(@(a) a.XLim(end), curAxes));
    curMaxY = max(arrayfun(@(a) a.YLim(end), curAxes));
    set(curAxes, 'XLim', [0, curMaxX], 'YLim', [0, curMaxY]);
    
    ylabel(curAxes(end), sprintf('%s axons', curAxonClass.tag));
end

xlabel(fig.Children(numTargetClasses), 'Synapses');

allAxes = flip(fig.Children);
arrayfun(@(ax, targetClass) title( ...
        ax, cat(1, targetClass, {'specific axons'}), ...
        'FontWeight', 'normal', 'FontSize', 10), ...
    allAxes(1:numTargetClasses), targetClasses);

annotation( ...
    fig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
	'String', {info.filename; info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');
