% This script simulates the synapse-seeded specificity analyses typically
% done manually by our experimentalists.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

minSynPre = 10;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, syn, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);

%% Prepare data
[conn, axonClasses] = ...
    connectEM.Connectome.prepareForSpecificityAnalysis( ...
        conn, axonClasses, 'minSynPre', minSynPre);
[classConn, targetClasses] = ...
    connectEM.Connectome.buildClassConnectome(conn);

% Build list of possible seed synapses
synT = connectEM.Connectome.buildSynapseTable(conn, syn);
synT.ontoTargetClass = conn.denMeta.targetClass(synT.postAggloId);
synT.isPriSpine = syn.synapses.type(synT.id) == 'PrimarySpine';

%% Build class with all axons
allAxonsClass = struct;
allAxonsClass.axonIds = find( ...
    conn.axonMeta.synCount >= minSynPre);
allAxonsClass.nullAxonIds = allAxonsClass.axonIds;
allAxonsClass.title = sprintf( ...
    'all axons with ≥ %d synapses (n = %d)', ...
    minSynPre, numel(allAxonsClass.axonIds));

axonClasses(end + 1) = allAxonsClass;

%% Seed synapses
colorT = table;
colorT.targetClass = targetClasses(:);
colors = get(groot, 'defaultAxesColorOrder');
colorT.color = colors(1:size(colorT, 1), :);

excSeedConfigs = { ...
    'ProximalDendrite', 'SmoothDendrite', 'ApicalDendrite'};
excSeedConfigs = cellfun(@(t) struct( ...
    'targetClass', t, ...
    'synIds', synT.id(synT.ontoTargetClass == t), ...
    'title', sprintf('Axons seeded at spine synapse onto %s', t), ...
    'color', colorT.color(colorT.targetClass == t, :)), ...
    reshape(excSeedConfigs, [], 1));

inhSeedConfigs = { ...
    'Somata', 'ProximalDendrite', 'SmoothDendrite', 'ApicalDendrite'};
inhSeedConfigs = cellfun(@(t) struct( ...
    'targetClass', t, ...
    'synIds', synT.id(synT.ontoTargetClass == t), ...
    'title', sprintf('Axons seeded at shaft synapse onto %s', t), ...
    'color', colorT.color(colorT.targetClass == t, :)), ...
    reshape(inhSeedConfigs, [], 1));

%% Axon populations
plotConfigs = axonClasses([1:4, end]);
[plotConfigs.seedConfigs] = deal(excSeedConfigs);
[plotConfigs([2, end]).seedConfigs] = deal(inhSeedConfigs);

%% Plot
for curIdx = 1:numel(plotConfigs)
    curConfig = plotConfigs(curIdx);
    withConfig(synT, classConn, targetClasses, info, false, curConfig);
    withConfig(synT, classConn, targetClasses, info, true, curConfig);
end

%% Functions
function withConfig(synT, classConn, targetClasses, info, weighted, config)
    fig = figure;
    fig.Color = 'white';
    fig.Position(3:4) = [600, 610];

    ax = axes(fig);
    ax.TickDir = 'out';
    hold(ax, 'on');
    
    nullSynCounts = sum(classConn(config.nullAxonIds, :), 2);
    nullSynProbs = sum(classConn(config.nullAxonIds, :), 1);
    nullSynProbs = nullSynProbs / sum(nullSynProbs);
    
    axonCounts = nan(size(config.seedConfigs));
    for curIdx = 1:numel(config.seedConfigs)
        curSeedConfig = config.seedConfigs(curIdx);
        curTargetClass = curSeedConfig.targetClass;
        
        curTargetClassId = find(targetClasses == curTargetClass);
        curNullSynProb = nullSynProbs(curTargetClassId);
        
       [obsConn, obsAxonIds, obsWeights] = forTargetClass( ...
            synT, classConn, targetClasses, curSeedConfig);
        
       [expConn, expAxonCounts] = ...
            connectEM.Specificity.calcExpectedDist( ...
                nullSynCounts, curNullSynProb, ...
                'distribution', 'binomial', ...
                'outputFormat', 'absolute');
            
        % Remove entries without seed synapse
        expAxonCounts(~expConn(:, 1)) = [];
        expConn(~expConn(:, 1), :) = [];
        
        % Remove seed synapse
        expSynCounts = expConn(:, 1);
        expConn = expConn - 1;
        
        % Restrict to axons of interest
        obsMask = ismember(obsAxonIds, config.axonIds);
        obsConn = obsConn(obsMask, :);
        obsWeights = obsWeights(obsMask);
        
        % Count axons (for legend)
        axonCounts(curIdx) = size(obsConn, 1);
        
        % Calculate and plot fractional synapses
        obsConn = obsConn ./ sum(obsConn, 2);
        expConn = expConn(:, 1) ./ expConn(:, 2);
        
        if weighted
            % Calculate weighted mean and standard deviation. Axons are
            % weighted by the probability of being reconstructed in sparse
            % synapse-seeded reconstructions.
            curWeights = obsWeights ./ sum(obsWeights);
            curExpWeights = expSynCounts .* expAxonCounts;
         
            % curExpWeights = curExpWeights .* expAxonCounts;
            curExpWeights = curExpWeights / sum(curExpWeights);
        else
            % Calculate unweighted mean and standard deviation.
            curWeights = ones(size(obsWeights)) / numel(obsWeights);
            curExpWeights = expAxonCounts / sum(expAxonCounts);
        end
        
        curMu = sum(curWeights .* obsConn, 1);
        curSigma = std(obsConn, curWeights, 1);
        
        curExpMu = sum(curExpWeights .* expConn);
        curExpSigma = std(expConn, curExpWeights);

        errorbar( ...
            1:numel(curMu), curMu, curSigma, ...
            'Color', curSeedConfig.color, ...
            'LineWidth', 1.25, ...
            'Marker', '.', ...
            'MarkerSize', 18);
        
        errorbar( ...
            curTargetClassId, curExpMu, curExpSigma, ...
            'Color', 'black', 'Marker', 'o', 'MarkerSize', 8);
    end
    
    plot( ...
        1:numel(nullSynProbs), nullSynProbs, ...
        'Color', 'black', ...
        'LineStyle', '--', ...
        'LineWidth', 2);

    xlim(ax, [0.5, (numel(targetClasses) + 0.5)]);
    xticks(ax, 1:numel(targetClasses));

    labels = arrayfun( ...
        @char, targetClasses, ...
        'UniformOutput', false);
    xticklabels(ax, labels);
    xtickangle(ax, 30);

    ylabel(ax, 'Fraction of synapses');
    ylim(ax, [0, 0.5]);
    
    legends = arrayfun( ...
        @(c, n) sprintf( ...
            '%s (n = %d)', c.title, n), ...
        config.seedConfigs, axonCounts, ...
        'UniformOutput', false);
    legends{end + 1} = 'Synapse fractions in null model';
    legends{end + 1} = 'Binomial null model';
    
    legend( ...
        flip(ax.Children([2, 1:2:end])), ...
        legends, 'Location', 'SouthOutside', 'Box', 'off');
    
    weightStr = {'unweighted', 'weighted'};
    weightStr = weightStr{1 + weighted};
    
    title(ax, ....
       {info.filename; info.git_repos{1}.hash; ...
        sprintf('%s (%s)', config.title, weightStr)}, ...
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
