% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
availFile = '/tmpscratch/amotta/l4/2018-02-02-surface-availability-connectome-axons-18-a/axon-avail-data.mat';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a.mat');

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

%% show how contact area translates into synapses
radiusUm = 1;
minSynCount = 10;

axonIds = find(conn.axonMeta.synCount >= minSynCount);
radiusIdx = find(avail.dists == (1E3 * radiusUm));

fig = figure();
for curClassIdx = 1:numel(targetClasses)
    curClassName = targetClasses(curClassIdx);
    
    curAvail = squeeze(availabilities( ...
        curClassIdx, radiusIdx, axonIds));
    curSpecs = squeeze(specificities( ...
        axonIds, curClassIdx));
    
    % plotting
    curAx = subplot(1, numel(targetClasses), curClassIdx);
    hold(curAx, 'on');
    
    % scatter plot
    scatter(curAx, curAvail, curSpecs, 'x');
    
    xlim(curAx, [0, 1]);
    ylim(curAx, [0, 1]);
    
    xlabel(curAx, ...
       {sprintf('A(%s)', curClassName); ...
        sprintf('at r_{pred} = %d µm', radiusUm)});
    ylabel(curAx, sprintf('S(%s)', curClassName));
    
    % diagonal
    plot(curAx, [0, 1], [0, 1], 'k--');
end

annotation( ...
    fig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'String', { ...
        'Target class availability versus specificity';
        info.git_repos{1}.hash});

fig.Position(3:4) = [1912, 290];

%% find AD specific axons
className = 'ApicalDendrite';
classIdx = find(targetClasses == className);
minSpecificity = 0.2;

axonIds = find( ...
    conn.axonMeta.synCount >= 10 ...
  & specificities(:, classIdx) >= minSpecificity);

%% plot
fig = figure();
ax = axes(fig);
hold(ax, 'on');

% color map
cmap = parula(101);
colormap(ax, cmap);
ax.Color = 'black';

maxSpec = max(specificities(axonIds, classIdx));
[~, sortIds] = sort(specificities(axonIds, classIdx), 'ascend');
axonIds = axonIds(sortIds);

for curAxonId = reshape(axonIds, 1, [])
    curAvails = shiftdim(availabilities(classIdx, :, curAxonId));
    
    curColor = specificities(curAxonId, classIdx);
    curColor = cmap(round(100 * curColor / maxSpec) + 1, :);
  
    plot(ax, avail.dists / 1E3, curAvails, 'Color', curColor);
end

xlabel(ax, 'r_{pred} (µm)');
ylabel(ax, sprintf('%s availability', className));
xlim(ax, [0, 60]);

cbar = colorbar('peer', ax);
cbar.Label.String = sprintf('%s specificity', className);
caxis([0, maxSpec]);

annotation( ...
    fig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'String', info.git_repos{1}.hash);

fig.Position(3:4) = [840, 630];

%% plot distribution over peak-predictability
% For each axon, we calculate
% * specificity / {min, max}(availability)
% * distribution over radii of peak predictability
classSpec = specificities(axonIds, classIdx);
classAvails = shiftdim(availabilities(classIdx, :, axonIds));

[classPeakAvail, classPeakIdx] = max(classAvails, [], 1);
classPeakPred = classSpec(:) ./ classPeakAvail(:);

fig = figure();
ax = axes(fig);

histogram(ax, classPeakPred);
ax.XLim = [0, 80];

xlabel(ax, { ...
    'specificity / max(availability)'; ...
    sprintf('for %s', className)});
ylabel(ax, {'Axons with'; sprintf( ...
    'S(%s) > %.1f', className, minSpecificity)});
title( ...
   {'Synapse predictability'; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);

%% plot availability vs. specificity
axonIds = find(conn.axonMeta.synCount >= 10);
peakAvails = squeeze(max(availabilities, [], 2))';

fig = figure();

for curClassIdx = 1:numel(targetClasses)
    curAx = subplot(1, numel(targetClasses), curClassIdx);
    scatter(curAx, ...
        specificities(axonIds, curClassIdx), ...
        peakAvails(axonIds, curClassIdx), 'x');
    
    xlim(curAx, [0, 1]);
    ylim(curAx, [0, 1]);
    
    xlabel(curAx, sprintf('S(%s)', targetClasses(curClassIdx)));
    ylabel(curAx, sprintf('A(%s)', targetClasses(curClassIdx)));
end