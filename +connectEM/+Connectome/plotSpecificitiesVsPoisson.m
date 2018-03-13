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
    fig.Position(3:4) = [1500, 300];
    
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
            specificities(:, classIdx), binEdges, ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2);

        xlabel(ax, sprintf('S(%s)', className));
        ax.XAxis.TickDirection = 'out';
        ax.XAxis.Limits = [0, 1];

        ax.YAxis.TickDirection = 'out';
        ax.YAxis.Limits(1) = 10 ^ (-0.1);
        ax.YAxis.Scale = 'log';

        axes{classIdx} = ax;
    end
    
    % Uncomment to show legend
    % leg = legend(ax, 'Expected (Poisson)', 'Observed');

    axes = horzcat(axes{:});
    yMax = max(arrayfun(@(a) a.YAxis.Limits(end), axes));
    for ax = axes; ax.YAxis.Limits(end) = yMax; end

    annotation(...
        'textbox', [0, 0.9, 1, 0.1], ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
        'String', ...
           {'Observed specificities vs. Poisson model'; ...
            axonClass.title; info.git_repos{1}.hash});
end