% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

% Partially split version of semi-automatically reconstructed axons,
% and dendrites with automatically and manually attached spines. This is
% the connectome used for the "L4 paper".
% connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

% Partially split version of automatically reconstructed axons,
% and dendrites with automatically attached spines. This is the
% automatically generated state prior to RoboEM.
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-04-partiallySplit-v1_dendrites-wholeCells-autoSpines-v1-classified-v2_SynapseAgglos-autoPreRobo-v1-classified.mat');

targetClasses = { ...
    'Somata', 'SO'; ...
    'ProximalDendrite', 'PD'; ...
    'ApicalDendrite', 'AD'; ...
    'SmoothDendrite', 'SD'; ...
    'AxonInitialSegment', 'AIS'; ...
    'OtherDendrite', 'OD'};

targetLabels = transpose(targetClasses(:, 2));
targetClasses = transpose(targetClasses(:, 1));

minSynPre = 10;
fdrThreshs = [0.025, 0.05, 0.1, 0.2, 0.3];

info = Util.runInfo();
Util.showRunInfo(info);

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, syn, axonClasses] = connectEM.Connectome.load(param, connFile);
[axonClasses.tag] = deal('Exc', 'Inh', 'TC', 'CC', 'Unclear');

% NOTE(amotta): Let's restrict to excitatory and inhibitory axons for now
axonClasses = axonClasses(1:2);

%% build class connectome for shaft synapses
% This is needed to simulate the effect of FP inhibitory synapses.
clear cur*;

curShaftConn = conn;
curShaftConn.connectome.synIdx = cellfun( ...
    @(ids) ids(syn.synapses.type(ids) == 'Shaft'), ...
    curShaftConn.connectome.synIdx, 'UniformOutput', false);

curShaftConn = ...
    connectEM.Connectome.prepareForSpecificityAnalysis( ...
        curShaftConn, axonClasses, 'minSynPre', minSynPre);
classConnShafts = ...
    connectEM.Connectome.buildClassConnectome( ...
        curShaftConn, 'targetClasses', targetClasses);
    
%% build complete class connectome
clear cur*;

[conn, axonClasses] = ...
    connectEM.Connectome.prepareForSpecificityAnalysis( ...
        conn, axonClasses, 'minSynPre', minSynPre);

classConn = ...
    connectEM.Connectome.buildClassConnectome( ...
        conn, 'targetClasses', targetClasses);

%% plot
clear cur*;

numRuns = numel(fdrThreshs);
evalT = cell(size(axonClasses));
barFracs = cell(size(axonClasses));

for curClassId = 1:numel(axonClasses)
    curAxonClass = axonClasses(curClassId);
    
    curEvalT = array2table( ...
        nan(numRuns, 1 + numel(targetClasses)), ...
        'VariableNames', [{'Overall'}, targetClasses]);
    curBarFracs = nan(numRuns, 2 * numel(targetClasses) - 1);
    
    for curRun = 1:numRuns
        curFdrThresh = fdrThreshs(curRun);
        curConn = classConn;
        
        % Only show plot for 20% false detection rate
        curShowPlot = (curFdrThresh == 0.2);
        
        curAxonClass.nullTargetClassProbs = ...
            connectEM.Specificity.calcFirstHitProbs( ...
                curConn(curAxonClass.nullAxonIds, :), ...
                'oneVersusRestBinomial');
        curAxonClass.specs =  ...
            connectEM.Specificity.calcForAxonClass( ...
                curConn, targetClasses, curAxonClass, ...
                'fdrThresh', curFdrThresh, ...
                'showPlot', curShowPlot, ...
                'info', info);
       [curA, curB, curC] = ...
            connectEM.Specificity.calcAxonFractions( ...
                curAxonClass, targetClasses);
        
        curEvalT{curRun, :} = [curA, curB];
        curBarFracs(curRun, :) = curC;
    end
    
    evalT{curClassId} = curEvalT;
    barFracs{curClassId} = curBarFracs;
end

%% quantitative evaluation
clear cur*;

for curIdx = 1:numel(axonClasses)
    curAxonClass = axonClasses(curIdx);
    
    curEvalT = evalT{curIdx};
    curEvalT.fdrThresh = fdrThreshs(:);
    curEvalT = circshift(curEvalT, 1, 2);
    
    fprintf('%s\n\n', curAxonClass.title);
    disp(curEvalT); fprintf('\n');
end

%% Plot fraction of axons specific
% Copy & paste from +connectEM/+Figure/fractionOfAxonsSpecific.m
clear cur*;
curPlotClasses = 1:2;

curFig = figure();
curAx = axes(curFig);
hold(curAx, 'on');

for curIdx = 1:numel(curPlotClasses)
    curId = curPlotClasses(curIdx);
    curFracs = evalT{curId}.Overall;
    
    plot(curAx, repelem(curIdx, numel(curFracs)), curFracs');
end

curLines = curAx.Children;
set(curLines, 'Color', 'black', 'Marker', '.', 'MarkerSize', 16);

curAx.XLim = [0.5, numel(curPlotClasses) + 0.5];
curAx.YLim = [0, 1];

curAx.XTick = 1:numel(curPlotClasses);
curAx.XTickLabel = {axonClasses(curPlotClasses).tag};
ylabel(curAx, {'Fraction of'; 'axons specific'});

connectEM.Figure.config(curFig, info);
curFig.Position(3:4) = [125, 160];

%% Plot distribution of specificities over target classes
clear cur*;

curBarWidth = 0.8;
curHistWidth = 0.09;
curPlotClasses = 1:2;

curColors = get(groot, 'defaultAxesColorOrder');
curColors = flip(curColors(1:numel(targetClasses), :), 1);

curMixed = (curColors(1:(end - 1), :) .^ 2 + curColors(2:end, :) .^ 2) / 2;
curMixed = cat(1, sqrt(curMixed), zeros(1, 3));

curColors = cat(1, transpose(curColors), transpose(curMixed));
curColors = transpose(reshape(curColors, 3, []));
curColors = num2cell(curColors(1:(end - 1), :), 2);

for curIdx = 1:numel(fdrThreshs)
    curFdrThresh = fdrThreshs(curIdx);
    
    curFig = figure();
    curFig.Color = 'white';
    curFig.Position(3:4) = [280, 280];

    curAx = axes(curFig); %#ok
    hold(curAx, 'on');
    curAx.YAxisLocation = 'right';

    curPlotData = cell2mat(cellfun( ...
        @(b) b(curIdx, :), barFracs(:), 'UniformOutput', false));
    curPlotData = flip(curPlotData(curPlotClasses, :), 2);

    curBars = bar(curAx, curPlotData, 'stacked', 'BarWidth', curBarWidth);
    set(curBars, {'FaceColor'}, curColors, 'EdgeColor', 'none');

    xlim(curAx, [0.5, numel(curPlotClasses) + 0.5]);
    xticklabels(curAx, {axonClasses(curPlotClasses).tag});
    xticks(curAx, 1:numel(curPlotClasses));

    curLeg = legend(curAx, ...
        flip(curBars(1:2:end)), targetLabels, ...
        'Location', 'EastOutside');
    title(curAx, sprintf('FDR threshold of %g', curFdrThresh));

    connectEM.Figure.config(curFig, info);
end

%% Calculate coinnervation matrix
clear cur*;
curFdThresh = 0.2;

for curAxonClass = axonClasses(2)
    curConn = classConn;
    
    curAxonClass.nullTargetClassProbs = ...
        connectEM.Specificity.calcFirstHitProbs( ...
            curConn(curAxonClass.nullAxonIds, :), ...
            'oneVersusRestBinomial');
    curAxonClass.specs =  ...
        connectEM.Specificity.calcForAxonClass( ...
            curConn, targetClasses, curAxonClass, ...
            'fdrThresh', curFdThresh, 'showPlot', false);
	
    if isfield(curAxonClass, 'synFalsePosRates') ...
            && not(isempty(curAxonClass.synFalsePosRates))
        % NOTE(amotta): If the user has specified synapse false positive
        % rates then let's corrected for the expectation. As discussed with
        % MH on 07.03.2019 we're using the unmodified class connectome to
        % defined the target-specific axons.
        curFpRates = curAxonClass.synFalsePosRates;
        
        curShaftConn = classConnShafts;
        curSpineConn = curConn - curShaftConn;
        
        curConn = ...
            curShaftConn .* (1 - curFpRates(1, :)) ...
          + curSpineConn .* (1 - curFpRates(2, :));
    end
	
    curSpecs = curAxonClass.specs;
    curSpecClasses = fieldnames(curSpecs);
    
   [~, curSpecClasses] = ismember(curSpecClasses, targetClasses);
    curSpecClasses = sort(curSpecClasses);
    
    curCoinMat = nan(numel(curSpecClasses), numel(targetClasses));
    curSpecClassRates = nan(1, numel(curSpecClasses));
    
    for curSpecIdx = 1:numel(curSpecClasses)
        curSpecClass = curSpecClasses(curSpecIdx);
        curAxonIds = curSpecs.(targetClasses{curSpecClass}).axonIds;
        
        curCoinVec = curConn(curAxonIds, :);
        curCoinVec(:, curSpecClass) = 0;
        curCoinVec = sum(curCoinVec, 1) / sum(curCoinVec(:));
        curCoinMat(curSpecIdx, :) = curCoinVec;
        
        curSpecClassRate = curConn(curAxonIds, :);
        curSpecClassRate = ...
            sum(curSpecClassRate(:, curSpecClass)) ...
         ./ sum(curSpecClassRate(:));
        curSpecClassRates(curSpecIdx) = curSpecClassRate;
    end
    
    % All axons
    curCoinVec = curConn(curAxonClass.axonIds, :);
    curCoinVec = sum(curCoinVec, 1) / sum(curCoinVec(:));
    curCoinMat = cat(1, curCoinMat, curCoinVec);
    
    curTargetLabels = targetLabels;
    assert(isequal(curTargetLabels{end}, 'OD'));
    curTargetLabels = curTargetLabels(1:(end - 1));
    
    plotIt(info, ...
        curTargetLabels, curAxonClass, curSpecClasses, ...
        curCoinMat, 'specClassRates', curSpecClassRates);
end

%% plotting
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
                && specClasses(curRow) == curCol
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
