% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');

targetClasses = { ...
    'Somata', 'ProximalDendrite', 'ApicalDendrite', ...
    'SmoothDendrite', 'AxonInitialSegment', 'OtherDendrite'};

minSynPre = 10;
info = Util.runInfo();

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, ~, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);

% Inhibitory whole cell → smooth dendrite
% Excitatory whole cell → proximal dendrite
wcMask = conn.denMeta.targetClass == 'WholeCell';
inMask = conn.denMeta.isInterneuron;

conn.denMeta.targetClass(wcMask &  inMask) = 'SmoothDendrite';
conn.denMeta.targetClass(wcMask & ~inMask) = 'ProximalDendrite';

classConnectome = ...
    connectEM.Connectome.buildClassConnectome( ...
        conn, 'targetClasses', targetClasses);
    
%% generate a class with all axons
allAxonClass = struct;
allAxonClass.axonIds = find( ...
    conn.axonMeta.synCount >= minSynPre);
allAxonClass.nullAxonIds = find( ...
    conn.axonMeta.synCount >= minSynPre);
allAxonClass.title = sprintf( ...
    'all axons with ≥ %d synapses (n = %d)', ...
    minSynPre, numel(allAxonClass.axonIds));

axonClasses(end + 1) = allAxonClass;

%% plot
for curIdx = 1:numel(axonClasses)
    plotAxonClass( ...
        info, conn.axonMeta, classConnectome, ...
        targetClasses, axonClasses(curIdx));
end

%% plotting
function plotAxonClass(info, axonMeta, classConn, targetClasses, axonClass)
    axonCount = numel(axonClass.axonIds);
    axonSpecs = classConn(axonClass.axonIds, :);
    axonSpecs = axonSpecs ./ sum(axonSpecs, 2);
    
    %% preparations
    axonNullProbs = connectEM.Specificity.calcChanceProbs( ...
        classConn, axonClass.axonIds, axonClass.nullAxonIds, ...
        'distribution', 'binomial');
    
    % calculate overall synapse probabilities
    targetClassSyns = sum(classConn(axonClass.nullAxonIds, :), 1);
    targetClassProbs = targetClassSyns / sum(targetClassSyns);
    
    %% plotting
    fig = figure;
    fig.Color = 'white';
    fig.Position(3:4) = [1850, 885];
    
    binEdges = linspace(0, 1, 21);
    axes = cell(size(targetClasses));
    pValAxes = cell(size(targetClasses));

    for classIdx = 1:numel(targetClasses)
        className = targetClasses{classIdx};
        classProb = targetClassProbs(classIdx);
        
        axonClassSpecs = axonSpecs(:, classIdx);
        axonClassNullProbs = axonNullProbs(:, classIdx);
        
        % Null hypothesis
       [nullSynFrac, nullAxonCount] = ...
            connectEM.Specificity.calcExpectedDist( ...
                axonMeta.synCount(axonClass.axonIds), ...
                classProb, 'distribution', 'binomial');
        
        nullBinId = discretize(nullSynFrac, binEdges);
        nullBinCount = accumarray(nullBinId, nullAxonCount);

        % Measured
        ax = subplot(3, numel(targetClasses), classIdx);
        axis(ax, 'square');
        hold(ax, 'on');
        
        histogram(ax, ...
            axonClassSpecs, ...
            'BinEdges', binEdges, ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2, ...
            'FaceAlpha', 1);
        histogram(ax, ...
            'BinEdges', binEdges, ...
            'BinCounts', nullBinCount, ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2, ...
            'FaceAlpha', 1);

        xlabel(ax, 'Synapse fraction');
        ax.XAxis.TickDirection = 'out';
        ax.XAxis.Limits = [0, 1];
        
        ylabel(ax, 'Axons');
        ax.YAxis.TickDirection = 'out';
        ax.YAxis.Limits(1) = 10 ^ (-0.1);
        ax.YAxis.Scale = 'log';
        
        title(ax, className, 'FontWeight', 'normal', 'FontSize', 10);
        axes{classIdx} = ax;
        
        %% p-values
        curBinEdges = linspace(-1E-3, 1 + 1E-3, numel(binEdges));
        
       [expChanceProbs, expChanceCounts] = ...
            connectEM.Specificity.calcExpectedChanceProbDist( ...
                axonMeta.synCount(axonClass.axonIds), classProb);
        curExpCounts = accumarray( ...
            discretize(expChanceProbs, curBinEdges), ...
            expChanceCounts / axonCount);
        curBinCounts = accumarray( ...
            discretize(axonClassNullProbs, curBinEdges), ...
            1 / axonCount);
            
        ax = subplot( ...
            3, numel(targetClasses), ...
            numel(targetClasses) + classIdx);
        axis(ax, 'square');
        hold(ax, 'on');
        
        histogram(ax, ...
            'BinCounts', curBinCounts, ...
            'BinEdges', curBinEdges, ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2);
        histogram(ax, ...
            'BinEdges', curBinEdges, ...
            'BinCounts', curExpCounts, ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2);
        
        ax.YScale = 'log';
        ax.XLim = curBinEdges([1, end]);
        
        pValAxes{classIdx} = ax;
        
        %% alternative visualization
        % Compare p-value distribution against expectation:
        % We'd expect there to be `theta` percent of axons with a p-value
        % below `theta`. If there are, however, significantly more axons
        % with a p-value below `theta`, something interesting is going on.
        curPVal = sort(axonClassNullProbs, 'ascend');
        curPVal = reshape(curPVal, 1, []);
        
        ax = subplot( ...
            3, numel(targetClasses), ...
            2 * numel(targetClasses) + classIdx);
        axis(ax, 'square');
        hold(ax, 'on');
        
       [curPVal, ~, curPAxonFrac] = unique(curPVal);
        curPAxonFrac = accumarray(curPAxonFrac, 1);
        curPAxonFrac = cumsum(curPAxonFrac) / sum(curPAxonFrac);
        
        curExpX = expChanceProbs;
        curExpY = cumsum(expChanceCounts);
        curExpY = curExpY / curExpY(end);
        
        curDiffs = interp1(curExpX, curExpY, curPVal);
        curDiffs = curPAxonFrac(:) - curDiffs(:);
        
        curThetaIdx = find(curDiffs(1:(end - 1)) < 0, 1);
       [curMaxDiff, curThetaIdx] = max(curDiffs(1:curThetaIdx));
        if curMaxDiff < 0; curThetaIdx = []; end
        
        plot(ax, curPVal, curPAxonFrac, 'LineWidth', 2);
        plot(ax, curExpX, curExpY, 'LineWidth', 2);
        
        if ~isempty(curThetaIdx)
            curThetaPVal = curPVal(curThetaIdx);
            plot(ax, ...
                repelem(curThetaPVal, 2), [0, 1], ...
                'Color', 'black', 'LineStyle', '--');
        end
        
        xlim(ax, [0, 1]);
        ylim(ax, [0, 1]);
        xlabel(ax, 'p-value');
        ylabel(ax, {'Fraction of axons'; 'with p < x'});
    end
    
    % Legend
    ax = axes{end};
    axPos = ax.Position;
    leg = legend(ax, ...
        'Observed', ...
        'Binomial model', ...
        'Location', 'East');
    leg.Box = 'off';
    
    % Fix positions
    ax.Position = axPos;
    leg.Position(1) = sum(axPos([1, 3])) + 0.005;

    axes = horzcat(axes{:});
    yMax = max(arrayfun(@(a) a.YAxis.Limits(end), axes));
    for ax = axes; ax.YAxis.Limits(end) = yMax; end
    
    annotation( ...
        'textbox', [0, 0.9, 1, 0.1], ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
        'String', { ...
            'Observed synapse fractions vs. null hypothesis'; ...
            axonClass.title; info.git_repos{1}.hash});
end
