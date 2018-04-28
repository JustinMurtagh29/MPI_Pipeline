% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');
availFile = '/tmpscratch/amotta/l4/2018-04-27-surface-availability-connectome-v5/axon-availability_v2.mat';

minSynPre = 10;
info = Util.runInfo();

%% loading data
avail = load(availFile);
conn = load(connFile);

%% build class connectome
% try to replicate cass connectome
connectome = conn.connectome;

% count synapses
connectome.synCount = cellfun(@numel, connectome.synIdx);
connectome.synIdx = [];

% add target class to connectome
[~, denMetaRow] = ismember(connectome.edges(:, 2), conn.denMeta.id);
connectome.targetClass = conn.denMeta.targetClass(denMetaRow);

% renumber classes alphabetically
conn.denMeta.targetClass = reordercats(conn.denMeta.targetClass);
targetClasses = unique(conn.denMeta.targetClass);

[~, connectome.targetClassId] = ...
    ismember(connectome.targetClass, targetClasses);

classConnectome = accumarray( ...
    cat(2, connectome.edges(:, 1), connectome.targetClassId), ...
    connectome.synCount, [numel(conn.axons), numel(targetClasses)]);
specificities = classConnectome ./ sum(classConnectome, 2);

%% build availabilities
[~, classIds] = ismember(targetClasses, avail.targetClasses);
availabilities = avail.axonAvail(classIds, :, :);
availabilities = availabilities ./ sum(availabilities, 1);

%% plot availability vs. distance
axonIds = find(conn.axonMeta.synCount >= minSynPre);
cmap = parula(101);

fig = figure();
fig.Color = 'white';

for classIdx = 1:numel(targetClasses)
    className = char(targetClasses(classIdx));
    
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
