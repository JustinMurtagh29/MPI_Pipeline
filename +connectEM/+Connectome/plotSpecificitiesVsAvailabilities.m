% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');
availFile = '/tmpscratch/amotta/l4/2018-04-27-surface-availability-connectome-v5/axon-availability_v2.mat';

minSynPre = 10;
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
        curColor = cmap(round(100 * curColor) + 1, :);

        plot(ax, avail.dists / 1E3, curAvails, 'Color', curColor);
    end

    xlabel(ax, {className; 'r_{pred} (µm)'});
    xlim(ax, [0, avail.dists(end)] / 1E3);
end

% Y axis
yMax = max(arrayfun(@(a) a.YLim(end), fig.Children));
[fig.Children.YLim] = deal([0, yMax]);
ylabel(fig.Children(end), 'Availability');

% X axis (as for panel 4c)
[fig.Children.XLim] = deal([0, 50]);

% colorbar
oldPos = ax.Position;
cbar = colorbar(ax);
cbar.Label.String = 'Specificity';
ax.Position = oldPos;

annotation( ...
    fig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'String', {info.filename; info.git_repos{1}.hash});
