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

%% loading data
conn = load(connFile);

%% sanity check
% try to replicate cass connectome
connectome = conn.connectome;

% remove dendrites with too few synapses
denMask = (conn.denMeta.synCount >= minSynPost);
connectome(~denMask(connectome.edges(:, 2)), :) = [];

% count synapses
connectome.synCount = cellfun(@numel, connectome.synIdx);
connectome.synIdx = [];

% add target class to connectome
[~, denMetaRow] = ismember(connectome.edges(:, 2), conn.denMeta.id);
connectome.targetClass = conn.denMeta.targetClass(denMetaRow);

% calculate # synapses per axon per target
targetClasses = { ...
    'Somata', 'WholeCell', 'ApicalDendrite', ...
    'SmoothDendrite', 'AxonInitialSegment', 'OtherDendrite'};
[~, connectome.targetClassId] = ismember( ...
    connectome.targetClass, targetClasses);

classConnectome = accumarray( ...
    cat(2, connectome.edges(:, 1), connectome.targetClassId), ...
    connectome.synCount, [numel(conn.axons), numel(targetClasses)]);
% den

%% specificity analysis
axonMask = (conn.axonMeta.synCount >= minSynPre);

specificities = classConnectome(axonMask, :);
specificities = specificities ./ sum(specificities, 2);

classCount = numel(conn.denClasses);

fig = figure;
for curIdx = 1:classCount
    curClassName = conn.denClasses{curIdx};
    
    ax = subplot(1, classCount, curIdx);
    histogram( ...
        ax, specificities(:, curIdx), ...
        linspace(0, 1, 51), 'EdgeColor', 'none');

    xlabel(ax, curClassName);
    ax.XAxis.TickDirection = 'out';
    ax.XAxis.Limits = [0, 1];
    
    if curIdx == 1
        ylabel(ax, 'Axons');
    else
        ax.YTickLabel = [];
    end
    
    ax.YAxis.TickDirection = 'out';
    ax.YAxis.Limits(1) = 10 ^ (-0.1);
    ax.YAxis.Scale = 'log';
end

annotation(...
    'textbox', [0, 0.9, 1, 0.1], ...
    'String', 'Target specificities', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')

yMax = max(arrayfun(@(a) a.YAxis.Limits(end), fig.Children));
for i = 1:numel(fig.Children); fig.Children(i).YAxis.Limits(end) = yMax; end

fig.Position(3:4) = [1099, 570];