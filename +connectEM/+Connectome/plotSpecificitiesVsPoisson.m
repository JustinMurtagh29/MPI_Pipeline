% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
param = struct;
param.saveFolder = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connName = 'connectome_axons_18_a_ax_spine_syn_clust';

targetClasses = { ...
    'Somata', 'WholeCell', 'ApicalDendrite', ...
    'SmoothDendrite', 'AxonInitialSegment', 'OtherDendrite'};

minSynPre = 10;
info = Util.runInfo();

%% loading data
conn = ...
    connectEM.Connectome.load(param, connName);
classConnectome = ...
    connectEM.Connectome.buildClassConnectome(conn, targetClasses);
axonClasses = ...
    connectEM.Connectome.buildAxonClasses(conn, 'minSynPre', minSynPre);

%% plot
for curIdx = 1:numel(axonClasses)
    plotAxonClass( ...
        info, conn.axonMeta, classConnectome, ...
        targetClasses, axonClasses(curIdx));
end

%% plotting
function plotAxonClass(info, axonMeta, classConn, targetClasses, axonClass)
    axonSpecs = classConn(axonClass.axonIds, :);
    axonSpecs = axonSpecs ./ sum(axonSpecs, 2);
    
    %% preparations
    % select axons for Poisson
    if isfield(axonClass, 'poissAxonIds') ...
            && ~isempty(axonClass.poissAxonIds)
        poissAxonIds = axonClass.poissAxonIds;
    else
        poissAxonIds = axonClass.axonIds;
    end
    
    axonPoissProbs = ...
        connectEM.Specificity.calcPoissonProbs( ...
            classConn, axonClass.axonIds, poissAxonIds);
    
    % calculate probabilities for Poisson
    targetClassSyns = sum(classConn(poissAxonIds, :), 1);
    targetClassProbs = targetClassSyns / sum(targetClassSyns);
    
    %% plotting
    fig = figure;
    fig.Color = 'white';
    fig.Position(3:4) = [1850, 400];
    
    binEdges = linspace(0, 1, 21);
    axes = cell(size(targetClasses));

    for classIdx = 1:numel(targetClasses)
        className = targetClasses{classIdx};
        classProb = targetClassProbs(classIdx);
        
        axonClassSpecs = axonSpecs(:, classIdx);
        axonClassPoissProbs = axonPoissProbs(:, classIdx);
        isSpecific = axonClassPoissProbs < 0.01;

        % Poisson
       [poissSynFrac, poissAxonCount] = ...
            connectEM.Specificity.calcExpectedDist( ...
                axonMeta.synCount(axonClass.axonIds), classProb, ...
                'distribution', 'poisson');
        
        poissBinId = discretize(poissSynFrac, binEdges);
        poissBinCount = accumarray(poissBinId, poissAxonCount);

        % Binomial
       [binoSynFrac, binoAxonCount] = ...
            connectEM.Specificity.calcExpectedDist( ...
                axonMeta.synCount(axonClass.axonIds), classProb, ...
                'distribution', 'binomial');
        
        binoBinId = discretize(binoSynFrac, binEdges);
        binoBinCount = accumarray(binoBinId, binoAxonCount);

        % Measured
        ax = subplot(1, numel(targetClasses), classIdx);
        axis(ax, 'square');
        hold(ax, 'on');
        
        histogram(ax, ...
            axonClassSpecs, ...
            'BinEdges', binEdges, ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2);
        histogram(ax, ...
            'BinEdges', binEdges, ...
            'BinCounts', poissBinCount, ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2);
        histogram(ax, ...
            'BinEdges', binEdges, ...
            'BinCounts', binoBinCount, ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2);

        xlabel(ax, className);
        ax.XAxis.TickDirection = 'out';
        ax.XAxis.Limits = [0, 1];

        ax.YAxis.TickDirection = 'out';
        ax.YAxis.Limits(1) = 10 ^ (-0.1);
        ax.YAxis.Scale = 'log';
        
        % Show number of overly specific axons. Two approaches are used:
        % 1. Find a synapse fraction threshold above which less than one
        %    axon is expected under the Poisson assumption. Then count the
        %    number of axons above this threshold.
        % 2. Count the number of axons which are above the Poisson
        %    distribution. This does not truly represent overly-specific
        %    axons.
        overBinId = 1 + find(poissBinCount > 1, 1, 'last');
        
        obsBinCount = histcounts(axonClassSpecs, binEdges);
        overThreshCount = sum(obsBinCount(overBinId:end));
        specificCount = sum(isSpecific);
        
        title(ax, { ...
            sprintf('%d with S ≥ %.1f', ...
                overThreshCount, binEdges(overBinId)); ...
            sprintf('%d with p ≤ 1 %%', specificCount)}, ...
            'FontWeight', 'normal', 'FontSize', 10);
        
        axes{classIdx} = ax;
    end
    
    ax = axes{end};
    axPos = ax.Position;
    leg = legend(ax, ...
        'Observed', ...
        'Poisson model', ...
        'Binomial model', ...
        'Location', 'East');
    
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
            'Synapse fractions vs. Poisson and binomial model'; ...
            axonClass.title; info.git_repos{1}.hash});
end