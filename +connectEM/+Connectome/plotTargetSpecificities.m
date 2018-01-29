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

%% spine fraction analysis
axonMask = (conn.axonMeta.synCount >= minSynPre);

% show spine synapse fraction for each axon
conn.axonMeta.spineSynFrac = ...
    conn.axonMeta.spineSynCount ...
    ./ conn.axonMeta.synCount;

fig = figure;
ax = axes(fig);

spineSynFrac = conn.axonMeta.spineSynFrac(axonMask, :);
histogram(ax, spineSynFrac, linspace(0, 1, 51));

xlim(ax, [0, 1]);
xlabel(ax, 'Fraction of synapses onto spines');
ylabel(ax, 'Axons');

ax.XAxis.TickDirection = 'out';
ax.YAxis.TickDirection = 'out';

fig.Position(3:4) = [570, 350];

%% build axon classes
axonClasses = { ...
    sprintf('Axons with spine fraction > %.1f', spineFracThresh)...
    axonMask & conn.axonMeta.spineSynFrac > spineFracThresh;
    sprintf('Axons with spine fraction < %.1f', spineFracThresh), ...
    axonMask & conn.axonMeta.spineSynFrac < spineFracThresh};

%% specificity analysis
preClassCount = size(axonClasses, 1);
postClassCount = numel(conn.denClasses);

fig = figure;
for curPreIdx = 1:preClassCount
    curPreName = axonClasses{curPreIdx, 1};
    curPreMask = axonClasses{curPreIdx, 2};
    
    specificities = classConnectome(curPreMask, :);
    specificities = specificities ./ sum(specificities, 2);

    for curPostIdx = 1:postClassCount
        curPostClassName = conn.denClasses{curPostIdx};
        
        curPlotIdx = curPostIdx + (curPreIdx - 1) * postClassCount;
        ax = subplot(preClassCount, postClassCount, curPlotIdx);
        
        histogram( ...
            ax, specificities(:, curPostIdx), ...
            linspace(0, 1, 51), 'EdgeColor', 'none');

        xlabel(ax, curPostClassName);
        ax.XAxis.TickDirection = 'out';
        ax.XAxis.Limits = [0, 1];

        if curPostIdx == 1
            ylabel(ax, curPreName);
        else
            ax.YTickLabel = [];
        end

        ax.YAxis.TickDirection = 'out';
        ax.YAxis.Limits(1) = 10 ^ (-0.1);
        ax.YAxis.Scale = 'log';
    end
end

annotation(...
    'textbox', [0, 0.9, 1, 0.1], ...
    'String', 'Target specificities', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')

yMax = max(arrayfun(@(a) a.YAxis.Limits(end), fig.Children));
for i = 1:numel(fig.Children); fig.Children(i).YAxis.Limits(end) = yMax; end

fig.Position(3:4) = [1099, 570];