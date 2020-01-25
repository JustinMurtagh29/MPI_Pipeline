% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
%   Benedikt Staffler <benedikt.staffler@brain.mpg.de>
%
% Differences to _v1:
%   - synapses: SynapseAgglos_v6_somaH.mat and the corresponding connectome

clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', ...
    'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v6-somaH.mat');

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

classConn = horzcat(classConn(:, targetIds), sum(classConn, 2));
classConn(:, end) = classConn(:, end) - sum(classConn(:, 1:(end - 1)), 2);

%% Calculate coinnervation matrix
for curAxonClass = axonClasses
    curSpecs = curAxonClass.specs;
    curSpecClasses = fieldnames(curSpecs);
    
   [~, curSpecClasses] = ismember(curSpecClasses, targetClasses);
    curSpecClasses = sort(curSpecClasses);
    
    curCoinMat = nan(1 + numel(curSpecClasses), 1 + numel(targetClasses));
    
    for curSpecIdx = 1:numel(curSpecClasses)
        curSpecClass = curSpecClasses(curSpecIdx);
        curAxonIds = curSpecs.(targetClasses{curSpecClass}).axonIds;
        
        curCoinVec = classConn(curAxonIds, :);
        curCoinVec(:, curSpecClass) = 0;
        curCoinVec = mean(curCoinVec ./ sum(curCoinVec, 2), 1);
        curCoinMat(curSpecIdx, :) = curCoinVec;
    end
    
    % All axons
    curCoinVec = classConn(curAxonClass.axonIds, :);
    curCoinVec = mean(curCoinVec ./ sum(curCoinVec, 2), 1);
    curCoinMat(end, :) = curCoinVec;
    
    plotIt(info, targetLabels, curAxonClass, curSpecClasses, curCoinMat);
end

%% Utilities
function plotIt(info, targetClasses, axonClass, specClasses, coinMat)
    maxDelta = 0.25;
    
    targetClasses{end + 1} = 'Other';
    specLabels = targetClasses(specClasses);
    specLabels{end + 1} = 'All';
    
    rows = numel(specLabels);
    cols = numel(targetClasses);
    frac = rows / cols;
    
    diagIds = arrayfun( ...
        @(idx, id) sub2ind(size(coinMat), idx, id), ...
        reshape(1:numel(specClasses), size(specClasses)), specClasses);
    
    deltaMat = coinMat - coinMat(end, :);
    deltaMat(diagIds) = 0;
    
    fig = figure();
    ax = axes(fig);
    
    imshow( ...
        deltaMat, [-maxDelta, +maxDelta], ...
        'Colormap', connectEM.Figure.blueWhiteRed(129), ...
        'Parent', ax);

    fig.Color = 'white';
    fig.Position(3:4) = 750 .* [1, frac];

    ax.Visible = 'on';
    ax.TickDir = 'out';
    ax.Box = 'off';

    ax.XAxisLocation = 'top';
    ax.XTick = 1:size(coinMat, 2);
    ax.XTickLabel = targetClasses;
    ax.XTickLabelRotation = 90;

    ax.YTick = 1:size(coinMat, 1);
    ax.YTickLabel = specLabels;
    ax.Position = [0.1, 0.01, 0.75, 0.75];
    
    for curIdx = 1:numel(coinMat)
       [curRow, curCol] = ind2sub(size(coinMat), curIdx);
        
        if curRow <= numel(specClasses) ...
                && specClasses(curRow) == curCol
            % Do not label "diagonals"
            continue;
        end

        curBoxSize = ax.Position(3:4) ./ [cols, rows];
        curOff = [curCol, numel(specLabels) - curRow + 1];
        curOff = ax.Position(1:2) + (curOff - 1) .* curBoxSize;
        
        curAnn = annotation( ...
            fig, 'textbox', [curOff, curBoxSize], ...
            'String', sprintf('%.1f %%', 100 * coinMat(curIdx)));
        curAnn.HorizontalAlignment = 'center';
        curAnn.VerticalAlignment = 'middle';
        curAnn.EdgeColor = 'none';
        curAnn.Color = 'black';
        curAnn.FontSize = 12;
    end
    
    cbar = colorbar('peer', ax);
    cbar.TickDirection = 'out';
    cbar.Position = [0.91, 0.1, 0.02, 0.8];
    cbar.Position([2, 4]) = ax.Position([2, 4]);
    cbar.TickLabels = arrayfun( ...
        @(f) sprintf('%+g %%', 100 * f), ...
        cbar.Ticks, 'UniformOutput', false);

    title(ax, ...
        {info.filename; info.git_repos{1}.hash; axonClass.title}, ...
        'FontWeight', 'normal', 'FontSize', 10);
end
