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
    specificities = classConn(axonClass.axonIds, :);
    specificities = specificities ./ sum(specificities, 2);
    
    %% preparations
    % select axons for Poisson
    if isfield(axonClass, 'poissAxonIds') ...
            && ~isempty(axonClass.poissAxonIds)
        poissAxonIds = axonClass.poissAxonIds;
    else
        poissAxonIds = axonClass.axonIds;
    end
    
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
    fig.Position(3:4) = [1500, 420];
    
    binEdges = linspace(0, 1, 11);
    axes = cell(size(targetClasses));

    for classIdx = 1:numel(targetClasses)
        className = targetClasses{classIdx};
        classProb = targetClassProbs(classIdx);

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
        obsBinCount = histcounts(specificities(:, classIdx), binEdges);

        % Measured
        ax = subplot(1, numel(targetClasses), classIdx);
        axis(ax, 'square');
        hold(ax, 'on');

        histogram(ax, ...
            'BinEdges', binEdges, ...
            'BinCounts', poissBinCount, ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2);

        histogram(ax, ...
            'BinEdges', binEdges, ...
            'BinCounts', obsBinCount, ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2);

        xlabel(ax, sprintf('S(%s)', className));
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
        overThreshCount = sum(obsBinCount(overBinId:end));
        overPoissCount = sum(max(obsBinCount(:) - poissBinCount(:), 0));
        
        overlyStr = { ...
            sprintf( ...
                '%d with S ≥ %.1f', ...
                overThreshCount, binEdges(overBinId)); ...
            sprintf( ...
                '%d above Poisson', round(overPoissCount))};
        
        annotation( ...
            fig, ...
            'textbox', ax.Position, ...
            'EdgeColor', 'none', ...
            'VerticalAlignment', 'bottom', ...
            'HorizontalAlignment', 'center', ...
            'String', overlyStr);
        
        axes{classIdx} = ax;
    end
    
    % Uncomment to show legend
    % leg = legend(ax, 'Expected (Poisson)', 'Observed');

    axes = horzcat(axes{:});
    yMax = max(arrayfun(@(a) a.YAxis.Limits(end), axes));
    for ax = axes; ax.YAxis.Limits(end) = yMax; end

    annotation( ...
        'textbox', [0, 0.9, 1, 0.1], ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
        'String', ...
           {'Observed specificities vs. Poisson model'; ...
            axonClass.title; info.git_repos{1}.hash});
end