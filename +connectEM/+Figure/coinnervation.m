% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a_ax_spine_syn_clust.mat');
synFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v3_ax_spine_clustered_classified.mat');

minSynCount = 10;
targetClasses = { ...
    'Somata', 'SO'; ...
    'ApicalDendrite', 'AD'; ...
    'SmoothDendrite', 'SD'; ...
    'AxonInitialSegment', 'AIS'};
targetLabels = targetClasses(:, 2);
targetClasses = targetClasses(:, 1);

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

conn = connectEM.Connectome.load(param, connFile, synFile);

%% Building class connectome
[classConn, targetIds] = ...
    connectEM.Connectome.buildClassConnectome(conn);
[~, targetIds] = ismember(targetClasses, targetIds);
assert(all(targetIds));

%% Calculate coinnervation matrix
axonIds = find(conn.axonMeta.synCount >= minSynCount);
coinMat = nan(numel(targetClasses));

for curIdx = 1:numel(targetClasses)
    curTargetId = targetIds(curIdx);
    
    curAxonIds = find(classConn(:, curTargetId));
    curAxonIds = intersect(curAxonIds, axonIds);
    
    % Calculate co-innervation
    % After removal of seed synapse
    curClassConn = classConn(curAxonIds, :);
    curClassConn(:, curTargetId) = curClassConn(:, curTargetId) - 1;
    
    curClassConn = curClassConn ./ sum(curClassConn, 2);
    curClassConn = curClassConn(:, targetIds);
    curClassConn = mean(curClassConn, 1);
    
    coinMat(curIdx, :) = curClassConn;
end

%% Plot matrix
% TODO(amotta): color
fig = figure();
ax = axes(fig);
imshow(10 .* coinMat, 'Parent', ax);

fig.Color = 'white';
fig.Position(3:4) = [660, 660];

ax.Visible = 'on';
ax.Box = 'off';
ax.TickDir = 'out';

ax.XAxisLocation = 'top';
ax.XTick = 1:numel(targetClasses);
ax.XTickLabel = targetLabels;
ax.XTickLabelRotation = 90;

ax.YTick = 1:numel(targetClasses);
ax.YTickLabel = targetLabels;

axis(ax, 'square');
ax.Position = [0.1, 0.1, 0.8, 0.8];

for curIdx = 1:numel(coinMat)
   [curRow, curCol] = ind2sub(size(coinMat), curIdx);
    
    curBoxSize = ax.Position(3:4) / numel(targetClasses);
    curOff = [curCol, numel(targetClasses) - curRow + 1];
    curOff = ax.Position(1:2) + (curOff - 1) .* curBoxSize;
    
    annotation( ...
        'textbox', [curOff, curBoxSize], ...
        'String', sprintf('%.1f %%', 100 * coinMat(curIdx)), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'EdgeColor', 'none', ...
        'Color', 'white');
end

title(ax, ...
    {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
