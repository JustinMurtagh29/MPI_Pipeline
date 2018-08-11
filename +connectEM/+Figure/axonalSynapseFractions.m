% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

targetClasses = { ...
    'Somata', 'SO'; ...
    'ProximalDendrite', 'PD'; ...
    'SmoothDendrite', 'SD'; ...
    'ApicalDendrite', 'AD'; ...
    'AxonInitialSegment', 'AIS'; ...
    'OtherDendrite', 'Other'};

% Show mean values for these axon classes
plotAxonClassIds = [1, 2];
plotAxonClassStyles = { ...
    {'Color', 'black', 'LineStyle', '-'}, ...
    {'Color', 'black', 'LineStyle', ':'}};

targetLabels = targetClasses(:, 2);
targetClasses = targetClasses(:, 1);

minSynPre = 10;
info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, syn, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);
[conn, axonClasses] = ...
    connectEM.Connectome.prepareForSpecificityAnalysis( ...
        conn, axonClasses, 'minSynPre', minSynPre);
[classConnectome, curSortIds] = ...
    connectEM.Connectome.buildClassConnectome(conn);
[~, curSortIds] = ismember(targetClasses, curSortIds);
classConnectome = classConnectome(:, curSortIds);

%% Calculate target class innervation probabilities for null model
clear cur*;

for curIdx = 1:numel(axonClasses)
    curNullProbs = classConnectome(axonClasses(curIdx).nullAxonIds, :);
    curNullProbs = connectEM.Specificity.calcFirstHitProbs(curNullProbs);
    axonClasses(curIdx).nullTargetClassProbs = curNullProbs;
end

%% Plot results
clear cur*;
curClassConn = ...
    classConnectome ...
 ./ sum(classConnectome, 2);

curAxonClassMeans = cell2mat(arrayfun(@(a)  ...
    mean(curClassConn(a.axonIds, :), 1), ...
    reshape(axonClasses(plotAxonClassIds), [], 1), ...
    'UniformOutput', false));

curClassConn = curClassConn( ...
    conn.axonMeta.synCount >= minSynPre, :);

boxGroups = 1:numel(targetClasses);
boxGroups = repmat(boxGroups, size(curClassConn, 1), 1);

% Add jitter
boxJitter = 0.5 * (rand(size(boxGroups)) - 0.5);
boxJitter = boxGroups + boxJitter;

fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [275, 300];

ax = axes(fig);
hold(ax, 'on');

% Fake plots for legend
colors = ax.ColorOrder(1:numel(targetClasses), :);
cellfun(@(c) plot(ax, nan, nan, '.', 'Color', c), num2cell(colors, 2));

scatter( ...
    ax, boxJitter(:), curClassConn(:), ...
    [], colors(boxGroups(:), :), '.');

for curTargetId = 1:numel(targetClasses)
    curX = curTargetId + [-0.45, +0.45];
    
    % Mean over all axons
    curMean = mean(curClassConn(:, curTargetId));
    curColor = colors(curTargetId, :);
    
    plot(ax, ...
        curX, repelem(curMean, 2), ...
        'Color', curColor, 'LineWidth', 2);
    
    % Mean over axon classes
    for curAxonIdx = 1:size(curAxonClassMeans, 1)
        curAxonClass = axonClasses(plotAxonClassIds(curAxonIdx));
        curNullProb = curAxonClass.nullTargetClassProbs(curTargetId);
        curStyle = plotAxonClassStyles{curAxonIdx};
        
        plot(ax, ...
            curX, repelem(curNullProb, 2), ...
            'LineWidth', 2, curStyle{:});
    end
end

ax.Box = 'off';
ax.TickDir = 'out';
ax.YScale = 'log';

yticklabels(ax, arrayfun( ...
    @(f) sprintf('%g', f), ...
    yticks(ax), 'UniformOutput', false));
ax.YAxis.MinorTick = 'off';
ylabel(ax, 'Synapse fraction');

xlim(ax, [0, numel(targetClasses)] + 0.5);
xticks(ax, 1:numel(targetClasses));
xticklabels(ax, targetLabels);

title(ax, ...
    {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
