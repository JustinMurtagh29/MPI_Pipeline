% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');
availFile = '/tmpscratch/amotta/l4/2018-04-27-surface-availability-connectome-v5/axon-availability_v2.mat';

minSynPre = 10;
maxAvail = 0.7;
maxSpec = 0.7;
maxDist = 50;

targetClasses = { ...
    'Somata', ...
    'ProximalDendrite', ...
    'SmoothDendrite', ...
    'ApicalDendrite', ...
    'AxonInitialSegment'};
    
info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

conn = connectEM.Connectome.load(param, connFile);
avail = load(availFile);

%% Prepare data
% Rename whole cells to proximal dendrites
% TODO(amotta): Use `connectEM.Connectome.prepareForSpecificityAnalysis`
% once availability calculation was updated accordingly!
wcMask = conn.denMeta.targetClass == 'WholeCell';
conn.denMeta.targetClass(wcMask) = 'ProximalDendrite';

wcMask = avail.targetClasses == 'WholeCell';
avail.targetClasses(wcMask) = 'ProximalDendrite';

[specificities, allTargetClasses] = ...
    connectEM.Connectome.buildClassConnectome(conn);

[~, classIds] = ismember(allTargetClasses, avail.targetClasses);
availabilities = avail.axonAvail(classIds, :, :);

% Calculat fraction
availabilities = availabilities ./ sum(availabilities, 1);
specificities = specificities ./ sum(specificities, 2);

% Restrict to target classes to plot
[~, classIds] = ismember(targetClasses, allTargetClasses);
availabilities = availabilities(classIds, :, :);
specificities = specificities(:, classIds);

%% plot availability vs. distance
axonIds = find(conn.axonMeta.synCount >= minSynPre);
cmap = parula(101);

fig = figure();
fig.Color = 'white';

for classIdx = 1:numel(targetClasses)
    className = targetClasses{classIdx};
    
    % Sort axons (most specific axons will be shown in front)
   [~, sortIds] = sort(specificities(axonIds, classIdx), 'ascend');
    axonIds = axonIds(sortIds);
    
    ax = subplot(1, numel(targetClasses), classIdx);
    
    ax.Color = 'white';
    ax.TickDir = 'out';
    
    colormap(ax, cmap);
    hold(ax, 'on');

    for curAxonId = reshape(axonIds, 1, [])
        curAvails = shiftdim(availabilities(classIdx, :, curAxonId));

        curColor = specificities(curAxonId, classIdx);
        curColor = (size(cmap, 1) - 1) * min(curColor, maxSpec) / maxSpec;
        curColor = cmap(round(curColor) + 1, :);

        plot(ax, avail.dists / 1E3, curAvails, 'Color', curColor);
    end
    
    xlabel(ax, 'r_{pred} (??m)');
    title(ax, className, 'FontWeight', 'normal', 'FontSize', 10);
end

% Y axis
[fig.Children.YLim] = deal([0, maxAvail]);
ylabel(fig.Children(end), 'Availability');

% X axis (as for panel 4c)
[fig.Children.XLim] = deal([0, maxDist]);

% Colorbar
oldPos = ax.Position;
cbar = colorbar(ax);

cbar.Ticks = (0:0.1:maxSpec) / maxSpec;
cbar.TickLabels = arrayfun( ...
    @num2str, 0:0.1:maxSpec, ...
    'UniformOutput', false);
cbar.Label.String = 'Synapse fraction';

ax.Position = oldPos;

annotation( ...
    fig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'String', {info.filename; info.git_repos{1}.hash});
