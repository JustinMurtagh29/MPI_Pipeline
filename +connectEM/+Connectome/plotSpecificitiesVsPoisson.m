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
conn = connectEM.Connectome.load(param, connName);

%% build class connectome
% try to replicate cass connectome
connectome = conn.connectome;

% count synapses
connectome.synCount = cellfun(@numel, connectome.synIdx);
connectome.synIdx = [];

% add target class to connectome
[~, denMetaRow] = ismember(connectome.edges(:, 2), conn.denMeta.id);
connectome.targetClass = conn.denMeta.targetClass(denMetaRow);

[~, connectome.targetClassId] = ismember( ...
    connectome.targetClass, targetClasses);

classConnectome = accumarray( ...
    cat(2, connectome.edges(:, 1), connectome.targetClassId), ...
    connectome.synCount, [numel(conn.axons), numel(targetClasses)]);

%% build axon classes
conn.axonMeta.spineSynFrac = ...
    conn.axonMeta.spineSynCount ...
    ./ conn.axonMeta.synCount;

axonClasses = struct;
axonClasses(1).axonIds = find( ...
    conn.axonMeta.synCount >= minSynPre ...
  & conn.axonMeta.spineSynFrac > 0.7);
axonClasses(1).title = sprintf( ...
   ['axons with ≥ %d synapses and ', ...
    'at least 70 %% onto spines (n = %d)'], ...
    minSynPre, numel(axonClasses(end).axonIds));

axonClasses(2).axonIds = find( ...
    conn.axonMeta.synCount >= minSynPre ...
  & conn.axonMeta.spineSynFrac < 0.3);
axonClasses(2).title = sprintf( ...
   ['axons with ≥ %d synapses and ', ...
    'at most 30 %% onto spines (n = %d)'], ...
    minSynPre, numel(axonClasses(end).axonIds));

axonClasses(3).axonIds = find( ...
    conn.axonMeta.isThalamocortical);
axonClasses(3).poissAxonIds = find( ...
    conn.axonMeta.synCount >= minSynPre ...
  & conn.axonMeta.spineSynFrac > 0.7);
axonClasses(3).title = sprintf( ...
    'thalamocortical axons (n = %d)', ...
    numel(axonClasses(end).axonIds));

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

    % restrict to axons with enough synapses
   [synCounts, ~, synCountAxons] = unique( ...
        axonMeta.synCount(axonClass.axonIds));
    synCountAxons = accumarray(synCountAxons, 1);
    
    %% plotting
    fig = figure;
    fig.Color = 'white';
    fig.Position(3:4) = [1850, 800];
    
    binEdges = linspace(0, 1, 21);
    axes = cell(size(targetClasses));
    pValAxes = cell(size(targetClasses));

    for classIdx = 1:numel(targetClasses)
        className = targetClasses{classIdx};
        classProb = targetClassProbs(classIdx);
        
        axonClassSpecs = axonSpecs(:, classIdx);
        axonClassPoissProbs = axonPoissProbs(:, classIdx);
        isSpecific = axonClassPoissProbs < 0.01;

        % Poissons
        poiss = table;
        poiss.prob = cell2mat(arrayfun(@(nSyn, nAxons) ...
            nAxons * poisspdf((0:nSyn)', nSyn * classProb), ...
            synCounts, synCountAxons, 'UniformOutput', false));
        poiss.spec = cell2mat(arrayfun( ...
            @(nSyn) (0:nSyn)' ./ nSyn, ...
            synCounts, 'UniformOutput', false));

        poiss.binId = discretize(poiss.spec, binEdges);
        poissBinCount = accumarray(poiss.binId, poiss.prob);

        % Measured
        ax = subplot(2, numel(targetClasses), classIdx);
        axis(ax, 'square');
        hold(ax, 'on');

        histogram(ax, ...
            axonClassSpecs(isSpecific), ...
            'BinEdges', binEdges, ...
            'EdgeColor', 'none', ...
            'FaceAlpha', 1);
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
        
        %% Distribution of p-values
        ax = subplot( ...
            2, numel(targetClasses), ...
            numel(targetClasses) + classIdx);
        axis(ax, 'square');
        hold(ax, 'on');
        
        % Plot cumulative distribution vs. diagonal
        curValX = sort(axonClassPoissProbs, 'ascend');
        curValX = reshape(curValX, 1, []);

        curValY = (1:numel(curValX)) ./ numel(curValX);
        curValY = curValY ./ curValX;
        
        curIdx = find(curValY > 1, 1);
        if mean(curValY(1:curIdx) > 1) < 0.5; curIdx = []; end
        
        if ~isempty(curIdx)
            curIdx = curIdx - 1 + find( ...
                curValY(curIdx:end) < 1, 1);
        end
        
        yyaxis(ax, 'right');
        hold(ax, 'on');
        plot(ax, curValX, curValY, 'LineWidth', 2);
        plot(ax, binEdges([1, end]), [1, 1], 'LineStyle', '--');
        ylim(ax, [0, 2]);
        
        if ~isempty(curIdx)
            plot(ax, curValX([curIdx, curIdx]), ax.YLim, 'k--');
            title(ax, ...
                sprintf('p = %.2f', curValX(curIdx)), ...
                'FontWeight', 'normal', 'FontSize', 10);
        end
        
        yyaxis(ax, 'left');
        histogram(ax, ...
            axonClassPoissProbs, binEdges, ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2);
        xlim(ax, binEdges([1, end]));
        ylim(ax, [0, size(axonSpecs, 1)]);
        
        xlabel(ax, 'p-value');
        ax.XAxis.TickDirection = 'out';
        
        pValAxes{classIdx} = ax;
    end
    
    % Uncomment to show legend
    % leg = legend(ax, 'Expected (Poisson)', 'Observed');

    axes = horzcat(axes{:});
    yMax = max(arrayfun(@(a) a.YAxis.Limits(end), axes));
    for ax = axes; ax.YAxis.Limits(end) = yMax; end
    
    pValAxes = horzcat(pValAxes{:});
    yMax = max(arrayfun(@(a) a.YAxis(1).Limits(end), pValAxes));
   [pValAxes.YLim] = deal([0, yMax]);

    annotation( ...
        'textbox', [0, 0.9, 1, 0.1], ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
        'String', ...
           {'Observed synapse fractions vs. Poisson model'; ...
            axonClass.title; info.git_repos{1}.hash});
end