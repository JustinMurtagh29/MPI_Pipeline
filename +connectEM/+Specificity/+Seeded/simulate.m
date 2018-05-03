% This script simulates the synapse-seeded specificity analyses typically
% done manually by our experimentalists.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');

minSynPre = 10;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, syn, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);

%% Prepare data
conn = ...
    connectEM.Connectome.prepareForSpecificityAnalysis(conn);
[classConn, targetClasses] = ...
    connectEM.Connectome.buildClassConnectome(conn);

% Build list of possible seed synapses
synT = connectEM.Connectome.buildSynapseTable(conn, syn);
synT.ontoTargetClass = conn.denMeta.targetClass(synT.postAggloId);

%% Seed synapses
colorT = table;
colorT.targetClass = targetClasses(:);
colors = get(groot, 'defaultAxesColorOrder');
colorT.color = colors(1:size(colorT, 1), :);

excSeedConfigs = { ...
    'ProximalDendrite', 'ApicalDendrite'};
excSeedConfigs = cellfun(@(t) struct( ...
    'synIds', synT.id(synT.isSpine & synT.ontoTargetClass == t), ...
    'title', sprintf('Axons seeded at spine synapse onto %s', t), ...
    'color', colorT.color(colorT.targetClass == t, :)), ...
    reshape(excSeedConfigs, [], 1));

inhSeedConfigs = { ...
    'Somata', 'ProximalDendrite', ...
    'ApicalDendrite', 'SmoothDendrite'};
inhSeedConfigs = cellfun(@(t) struct( ...
    'synIds', synT.id(~synT.isSpine & synT.ontoTargetClass == t), ...
    'title', sprintf('Axons seeded at shaft synapse onto %s', t), ...
    'color', colorT.color(colorT.targetClass == t, :)), ...
    reshape(inhSeedConfigs, [], 1));

%% Axon populations
plotConfigs = axonClasses([4, 3, 2]);
plotConfigs = rmfield(plotConfigs, 'nullAxonIds');

[plotConfigs(1:2).seedConfigs] = deal(excSeedConfigs);
[plotConfigs(3).seedConfigs] = deal(inhSeedConfigs);

%% Plot
for curIdx = 1:numel(plotConfigs)
    curConfig = plotConfigs(curIdx);
    withConfig(synT, classConn, targetClasses, info, curConfig);
end

%% Functions
function withConfig(synT, classConn, targetClasses, info, config)
    fig = figure;
    fig.Color = 'white';
    fig.Position(3:4) = [600, 610];

    ax = axes(fig);
    ax.TickDir = 'out';
    axis(ax, 'square');
    hold(ax, 'on');

    axonCounts = nan(size(config.seedConfigs));
    for curIdx = 1:numel(config.seedConfigs)
        curSeedConfig = config.seedConfigs(curIdx);
       [obsConn, obsAxonIds, obsWeights] = ...
            forTargetClass(synT, classConn, targetClasses, curSeedConfig);
        
        % Restrict to axons of interest
        obsMask = ismember(obsAxonIds, config.axonIds);
        obsConn = obsConn(obsMask, :);
        obsWeights = obsWeights(obsMask);
        
        % Count axons (for legend)
        axonCounts(curIdx) = size(obsConn, 1);
        
        % Calculate and plot fractional synapses
        obsConn = obsConn ./ sum(obsConn, 2);
        
        % Calculate (weighted) mean and standard deviation. Axons are
        % weighted by the probability of being reconstructed in sparse
        % synapse-seeded reconstructions.
        curWeights = obsWeights ./ sum(obsWeights);
        curMu = sum(curWeights .* obsConn, 1);
        curSigma = std(obsConn, curWeights, 1);

        errorbar( ...
            1:numel(curMu), curMu, curSigma, ...
            'Color', curSeedConfig.color, ...
            'LineWidth', 1.25, ...
            'Marker', '.', ...
            'MarkerSize', 18);
    end

    xlim(ax, [0.5, (numel(targetClasses) + 0.5)]);
    xticks(ax, 1:numel(targetClasses));

    labels = arrayfun( ...
        @char, targetClasses, ...
        'UniformOutput', false);
    xticklabels(ax, labels);
    xtickangle(ax, 30);

    ylabel(ax, 'Fraction of synapses');
    ylim(ax, [0, 1]);
    
    legends = arrayfun( ...
        @(c, n) sprintf( ...
            '%s (n = %d)', c.title, n), ...
        config.seedConfigs, axonCounts, ...
        'UniformOutput', false);
    legend(legends, 'Location', 'North', 'Box', 'off');
    
    title(ax, ....
       {info.filename; info.git_repos{1}.hash; config.title}, ...
        'FontWeight', 'normal', 'FontSize', 10);
end

function [obsConn, axonIds, axonSeedCount] = ...
        forTargetClass(seedSynT, classConn, targetClasses, seedConfig)
    % Determine seed synapses
    seedSynMask = ismember(seedSynT.id, seedConfig.synIds);
    seedSynT = seedSynT(seedSynMask, :);
    
    % Identified seeded axons
   [axonIds, uniRows] = unique(seedSynT.preAggloId);
    seedSynT = seedSynT(uniRows, :);
    
    % Trace axons
    obsConn = classConn(axonIds, :);
    
    % Remove seed synapses
   [~, slotIds] = ismember( ...
        seedSynT.ontoTargetClass, targetClasses);
    slotIds = arrayfun( ...
        @(r, c) sub2ind(size(obsConn), r, c), ...
        transpose(1:size(seedSynT, 1)), slotIds);
    obsConn(slotIds) = obsConn(slotIds) - 1;
    
    % Determine the number of synapses each axon made onto the seed class.
    % This information is needed to simulate the seeding probability in
    % sparse reconstructions.
    axonSeedCount = obsConn(slotIds) + 1;
end
