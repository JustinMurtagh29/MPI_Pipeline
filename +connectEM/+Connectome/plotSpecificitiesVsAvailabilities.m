% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
availFile = '/tmpscratch/amotta/l4/2018-02-02-surface-availability-connectome-axons-18-a/axon-avail-data.mat';

connFile = fullfile( ...
    rootDir, 'connectomeState', ...
    'connectome_axons_18_a_with_den_meta.mat');

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

%% build availabilitie
[~, classIds] = ismember(targetClasses, avail.targetClasses);
availabilities = avail.axonAvail(classIds, :, :);
availabilities = availabilities ./ sum(availabilities, 1);

%% find AD specific axons
className = 'ApicalDendrite';
classIdx = find(targetClasses == className);
minSpecificity = 0.25;

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

for curAxonId = reshape(axonIds, 1, [])
    curAvails = shiftdim(availabilities(classIdx, :, curAxonId));
    
    curColor = specificities(curAxonId, classIdx);
    curColor = cmap(round(100 * curColor) + 1, :);
  
    plot(ax, avail.dists / 1E3, curAvails, 'Color', curColor);
end

xlabel(ax, 'r_{pred}');
xlim(ax, [0, avail.dists(end)] / 1E3);
ylabel(ax, sprintf('%s availability', className));

cbar = colorbar('peer', ax);
cbar.Label.String = sprintf('%s specificity', className);

annotation( ...
    fig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'String', info.git_repos{1}.hash);