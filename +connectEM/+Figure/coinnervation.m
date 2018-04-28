% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');

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

[conn, ~, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);

%% Building class connectome
[classConn, targetIds] = ...
    connectEM.Connectome.buildClassConnectome(conn);
[~, targetIds] = ismember(targetClasses, targetIds);
assert(all(targetIds));

%% Calculate coinnervation matrix
axonIds = find(conn.axonMeta.synCount >= minSynCount);
coinMat = nan(numel(targetClasses));

for curIdx = 1:numel(targetClasses)
    curTargetClass = targetClasses{curIdx};
    curTargetId = targetIds(curIdx);
    
    % Collect specific axons
    curAxonIds = zeros(0, 1);
    for curAxonClass = reshape(axonClasses, 1, [])
        try
            curSpecs = curAxonClass.specs;
            curSpecs = curSpecs.(curTargetClass);
            curAxonIds = [curAxonIds; curSpecs.axonIds];
        catch
            % Yes, I'm abusing exceptions for control flow.
            % Whatcha gonna do 'bout it, ha? - amotta
        end
    end
    
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
imshow(coinMat, 'Parent', ax, 'Colormap', parula);

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
    
    curAnn = annotation( ...
        'textbox', [curOff, curBoxSize], ...
        'String', sprintf('%.1f %%', 100 * coinMat(curIdx)));
    curAnn.HorizontalAlignment = 'center';
	curAnn.VerticalAlignment = 'middle';
	curAnn.EdgeColor = 'none';
	curAnn.Color = 'white';
    curAnn.FontSize = 12;
end

cbar = colorbar('peer', ax);
cbar.TickDirection = 'out';
cbar.Position = [0.925, 0.1, 0.02, 0.8];

title(ax, ...
    {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
