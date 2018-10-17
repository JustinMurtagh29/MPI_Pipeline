% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

targetClasses = { ...
    'Somata', 'SO'; ...
    'ProximalDendrite', 'PD'; ...
    'ApicalDendrite', 'AD'; ...
    'SmoothDendrite', 'SD'; ...
    'AxonInitialSegment', 'AIS'};
targetLabels = targetClasses(:, 2);
targetClasses = targetClasses(:, 1);

minSynPre = 10;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, ~, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);
[conn, axonClasses] = ...
    connectEM.Connectome.prepareForSpecificityAnalysis( ...
        conn, axonClasses, 'minSynPre', minSynPre);

%% Find specific exc. and inh. axons
axonClasses = axonClasses(1:2);

axonClasses = ...
    connectEM.Connectome.buildAxonSpecificityClasses(conn, axonClasses);
    
%% Build class connectome
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
    curSpecClassRates = nan(1, numel(curSpecClasses));
    
    for curSpecIdx = 1:numel(curSpecClasses)
        curSpecClass = curSpecClasses(curSpecIdx);
        curAxonIds = curSpecs.(targetClasses{curSpecClass}).axonIds;
        
        curCoinVec = classConn(curAxonIds, :);
        curCoinVec(:, curSpecClass) = 0;
        curCoinVec = sum(curCoinVec, 1) / sum(curCoinVec(:));
        curCoinMat(curSpecIdx, :) = curCoinVec;
        
        curSpecClassRate = classConn(curAxonIds, :);
        curSpecClassRate = ...
            sum(curSpecClassRate(:, curSpecClass)) ...
         ./ sum(curSpecClassRate(:));
        curSpecClassRates(curSpecIdx) = curSpecClassRate;
    end
    
    % All axons
    curCoinVec = classConn(curAxonClass.axonIds, :);
    curCoinVec = sum(curCoinVec, 1) / sum(curCoinVec(:));
    curCoinMat(end, :) = curCoinVec;
    
    plotIt(info, ...
        targetLabels, curAxonClass, curSpecClasses, curCoinMat, ...
        'specClassRates', curSpecClassRates);
end

%% Utilities
function plotIt( ...
        info, targetClasses, axonClass, ...
        specClasses, coinMat, varargin)
    opt = struct;
    opt.maxDelta = 0.20;
    opt.specClassRates = [];
    opt = Util.modifyStruct(opt, varargin{:});
    
    specLabels = targetClasses(specClasses);
    specLabels{end + 1} = 'All';
    
    % Don't show the synapses, which do not fall in one of the target
    % classes listed in `targetClasses`.
    coinMat(:, end) = [];
    
    rows = numel(specLabels);
    cols = numel(targetClasses);
    frac = rows / cols;
    
    % Sanity checks
    assert(isequal(size(coinMat), [rows, cols]));
    
    diagIds = arrayfun( ...
        @(idx, id) sub2ind(size(coinMat), idx, id), ...
        reshape(1:numel(specClasses), size(specClasses)), specClasses);
    
    deltaMat = coinMat - coinMat(end, :);
    deltaMat(diagIds) = 0;
    
    if ~isempty(opt.specClassRates)
        coinMat(diagIds) = opt.specClassRates;
    end
    
    fig = figure();
    ax = axes(fig);
    
    imshow( ...
        deltaMat, [-opt.maxDelta, +opt.maxDelta], ...
        'Colormap', buildColormap(129), ...
        'Parent', ax);

    fig.Color = 'white';
    fig.Position(3:4) = 350 .* [1, frac];

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
        curEdgeColor = 'none';
        
        if curRow <= numel(specClasses) ...
                && specClasses(curRow) == curCol ...
                
            if isempty(opt.specClassRates); continue; end
            curEdgeColor = 'black';
        end

        curBoxSize = ax.Position(3:4) ./ [cols, rows];
        curOff = [curCol, numel(specLabels) - curRow + 1];
        curOff = ax.Position(1:2) + (curOff - 1) .* curBoxSize;
        
        curAnn = annotation( ...
            fig, 'textbox', [curOff, curBoxSize], ...
            'String', sprintf('%.2g', 100 * coinMat(curIdx)));
        curAnn.HorizontalAlignment = 'center';
        curAnn.VerticalAlignment = 'middle';
        curAnn.EdgeColor = curEdgeColor;
        curAnn.Color = 'black';
        curAnn.FontSize = 12;
        curAnn.LineWidth = 2;
    end
    
    cbar = colorbar('peer', ax);
    cbar.TickLabels = arrayfun( ...
        @(f) sprintf('%+g %%', 100 * f), ...
        cbar.Ticks, 'UniformOutput', false);
    cbar.TickDirection = 'out';
    cbar.Position = [0.91, 0.1, 0.02, 0.8];
    cbar.Position([2, 4]) = ax.Position([2, 4]);

    title(ax, ...
        {info.filename; info.git_repos{1}.hash; axonClass.title}, ...
        'FontWeight', 'normal', 'FontSize', 10);
end

function cmap = buildColormap(n)
    c = 1 + (n - 1) / 2;
    alpha = linspace(0, 1, c);
    alpha = transpose(alpha);
    
    cmap = zeros(n, 3);
    cmap(1:c, :) = alpha .* [1, 1, 1];
    cmap(c:n, :) = sqrt(( ...
        alpha .* [0.301, 0.745, 0.933] .^ 2 ...
      + (1 - alpha) .* [1, 1, 1] .^ 2));
end
