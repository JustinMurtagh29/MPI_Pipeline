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
synT.type = syn.synapses.type(synT.id);

%% Build class with all axons
clear cur*;

curAllAxonClass = struct;
curAllAxonClass.axonIds = find( ...
    conn.axonMeta.synCount >= minSynPre);
curAllAxonClass.nullAxonIds = curAllAxonClass.axonIds;
curAllAxonClass.title = sprintf( ...
    'all axons with ≥ %d synapses (n = %d)', ...
    minSynPre, numel(curAllAxonClass.axonIds));

axonClasses(end + 1) = curAllAxonClass;

%% Seed synapses
clear cur*;

curColorT = table;
curColorT.targetClass = targetClasses(:);
curColors = get(groot, 'defaultAxesColorOrder');
curColorT.color = curColors(1:size(curColorT, 1), :);

curInhSeedTargets = { ...
    'Somata', 'ProximalDendrite', 'SmoothDendrite', 'ApicalDendrite'};
curSeedConfigs = cellfun( ...
    @(t) struct( ...
        'targetClass', t, ...
        'synIds', synT.id(synT.ontoTargetClass == t), ...
        'title', sprintf('Axons seeded at synapse onto %s', t), ...
        'color', curColorT.color(curColorT.targetClass == t, :)), ...
    reshape(curInhSeedTargets, [], 1));

% inhibitory axons only
plotConfigs = axonClasses(2);
plotConfigs.seedConfigs = curSeedConfigs;

%% Plot
clear cur*;

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
    
    axonCounts = nan(size(config.seedConfigs));
    for curIdx = 1:numel(config.seedConfigs)
        curSeedConfig = config.seedConfigs(curIdx);
      	% curTargetClass = curSeedConfig.targetClass;
        
       [obsConn, obsAxonIds, obsWeights] = forTargetClass( ...
            synT, classConn, targetClasses, curSeedConfig);
        
        % Restrict to axons of interest
        obsMask = ismember(obsAxonIds, config.axonIds);
        obsConn = obsConn(obsMask, :);
        obsWeights = obsWeights(obsMask);
        
        % Count axons (for legend)
        axonCounts(curIdx) = size(obsConn, 1);
        
        % Calculate and plot fractional synapses
        obsConn = obsConn ./ sum(obsConn, 2);
        
        if weighted
            % Calculate weighted mean and standard deviation. Axons are
            % weighted by the probability of being reconstructed in sparse
            % synapse-seeded reconstructions.
            curWeights = obsWeights;
        else
            % Calculate unweighted mean and standard deviation.
            curWeights = ones(size(obsWeights));
        end
        
        curPerc = repelem(obsConn, curWeights, 1);
        curPerc = prctile(curPerc, [25, 50, 75], 1);
        
        curMed = curPerc(2, :);
        curNeg = curMed - curPerc(1, :);
        curPos = curPerc(3, :) - curMed;
        % curMu = sum(curWeights .* obsConn, 1);
        % curSigma = std(obsConn, curWeights, 1);

        errorbar( ...
            1:size(obsConn, 2), curMed, curNeg, curPos, ...
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
    ylim(ax, [0, 0.75]);
    
    legends = arrayfun( ...
        @(c, n) sprintf( ...
            '%s (n = %d)', c.title, n), ...
        config.seedConfigs, axonCounts, ...
        'UniformOutput', false);
    
    legend(ax, legends, 'Location', 'SouthOutside', 'Box', 'off');
    
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
