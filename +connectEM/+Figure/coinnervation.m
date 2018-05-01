% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');

targetClasses = { ...
    'Somata', 'SO'; ...
    'ProximalDendrite', 'PD'; ...
    'ApicalDendrite', 'AD'; ...
    'SmoothDendrite', 'SD'; ...
    'AxonInitialSegment', 'AIS'};
targetLabels = targetClasses(:, 2);
targetClasses = targetClasses(:, 1);

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, ~, axonClasses] = connectEM.Connectome.load(param, connFile);

%% Prepare connectome
conn = ...
    connectEM.Connectome.prepareForSpecificityAnalysis(conn);
axonClasses = ...
    connectEM.Connectome.buildAxonSpecificityClasses(conn, axonClasses);

[classConn, targetIds] = ...
    connectEM.Connectome.buildClassConnectome(conn);
[~, targetIds] = ismember(targetClasses, targetIds);
assert(all(targetIds));

synCount = sum(classConn, 2);
classConn = classConn(:, targetIds);
classConn(:, end + 1) = synCount - sum(classConn, 2);
classConn = classConn ./ sum(classConn, 2);

%% Calculate coinnervation matrix
for curAxonClass = axonClasses
    curSpecs = curAxonClass.specs;
    curSpecClasses = fieldnames(curSpecs);
    
   [~, curSpecClasses] = ismember(curSpecClasses, targetClasses);
    curSpecLabels = targetLabels(sort(curSpecClasses));
    curSpecClasses = targetClasses(sort(curSpecClasses));
    
    curCoinMat = nan(1 + numel(curSpecClasses), 1 + numel(targetClasses));
    
    for curSpecIdx = 1:numel(curSpecClasses)
        curSpecClass = curSpecClasses{curSpecIdx};
        curAxonIds = curSpecs.(curSpecClass).axonIds;
        
        curCoinVec = mean(classConn(curAxonIds, :), 1);
        curCoinMat(curSpecIdx, :) = curCoinVec;
    end
    
    % Other axons
    curAxonIds = structfun(@(curSpec) {curSpec.axonIds}, curSpecs);
    curAxonIds = setdiff(curAxonClass.axonIds, cell2mat(curAxonIds));
    curCoinMat(end, :) = mean(classConn(curAxonIds, :), 1);
    
    plotIt(info, targetLabels, curAxonClass, curSpecLabels, curCoinMat);
end

%% Utilities
function plotIt(info, targetClasses, axonClass, specClasses, coinMat)
    targetClasses{end + 1} = 'Other';
    specClasses{end + 1} = 'None';
    
    rows = numel(specClasses);
    cols = numel(targetClasses);
    frac = rows / cols;
    
    fig = figure();
    ax = axes(fig);
    imshow(coinMat, 'Parent', ax, 'Colormap', parula);

    fig.Color = 'white';
    fig.Position(3:4) = 750 .* [1, frac];

    ax.Visible = 'on';
    ax.Box = 'off';
    ax.TickDir = 'out';

    ax.XAxisLocation = 'top';
    ax.XTick = 1:size(coinMat, 2);
    ax.XTickLabel = targetClasses;
    ax.XTickLabelRotation = 90;

    ax.YTick = 1:size(coinMat, 1);
    ax.YTickLabel = specClasses;
    ax.Position = [0.1, 0.01, 0.75, 0.75];
    
    for curIdx = 1:numel(coinMat)
       [curRow, curCol] = ind2sub(size(coinMat), curIdx);

        curBoxSize = ax.Position(3:4) ./ [cols, rows];
        curOff = [curCol, numel(specClasses) - curRow + 1];
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
    cbar.Position([2, 4]) = ax.Position([2, 4]);

    title(ax, ...
        {info.filename; info.git_repos{1}.hash; axonClass.title}, ...
        'FontWeight', 'normal', 'FontSize', 10);
end
