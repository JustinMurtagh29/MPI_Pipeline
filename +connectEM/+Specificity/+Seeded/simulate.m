% This script simulates the synapse-seeded specificity analyses typically
% done manually by our experimentalists.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
param = struct;
param.saveFolder = '/gaba/u/mberning/results/pipeline/20170217_ROI';

connName = 'connectome_axons_18_a_ax_spine_syn_clust';
synFile = fullfile(param.saveFolder, 'connectomeState', 'SynapseAgglos_v3_ax_spine_clustered.mat');

minSynPre = 10;

info = Util.runInfo();

%% Loading data
syn = load(synFile);
conn = connectEM.Connectome.load(param, connName);

%% Process connectome
[classConn, targetClasses] = ...
    connectEM.Connectome.buildClassConnectome(conn);
axonClasses = ...
    connectEM.Connectome.buildAxonClasses(conn, 'minSynPre', minSynPre);

% Build list of possible seed synapses
synT = connectEM.Connectome.buildSynapseTable(conn, syn);
synT.ontoTargetClass = conn.denMeta.targetClass(synT.postAggloId);

%% Seed synapses
seedConfigs = struct;
seedConfigs(1).synIds = synT.id( ...
    ~synT.isSpine & synT.ontoTargetClass == 'ApicalDendrite');
seedConfigs(1).title = 'Axons seeded at shaft synapse onto apical dendrite';

seedConfigs(2).synIds = synT.id( ...
    ~synT.isSpine & synT.ontoTargetClass == 'Somata');
seedConfigs(2).title = 'Axons seeded at soma synapse';

%% Axon populations
plotConfigs = struct;
plotConfigs(1).axonIds = find( ...
    conn.axonMeta.synCount >= minSynPre);
plotConfigs(1).seedConfigs = seedConfigs;
plotConfigs(1).title = sprintf( ...
    'Axons with ≥ %d synapses', minSynPre);

plotConfigs(2).axonIds = union( ...
    axonClasses(1).axonIds, axonClasses(2).axonIds);
plotConfigs(2).seedConfigs = seedConfigs;
plotConfigs(2).title = sprintf( ...
    'Strongly exc. or inh. axons with ≥ %d synapses', minSynPre);

plotConfigs(3).axonIds = axonClasses(2).axonIds;
plotConfigs(3).seedConfigs = seedConfigs;
plotConfigs(3).title = sprintf( ...
    'Strongly inh. axons with ≥ %d synapses', minSynPre);

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
    
    titleString = sprintf( ...
        '%s (n = %d)', config.title, ...
        numel(config.axonIds));
    title(ax, ....
       {info.filename; info.git_repos{1}.hash; titleString}, ...
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