% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');
ctrlConnFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-linearized_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, syn, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);
[ctrlConn, ctrlSyn] = ...
    connectEM.Connectome.load(param, ctrlConnFile);

[~, synToSynFile] = fileparts(connFile);
synToSynFile = sprintf('%s_synToSynDists.mat', synToSynFile);
synToSynFile = fullfile(fileparts(connFile), synToSynFile);
synToSyn = load(synToSynFile);

[~, ctrlSynToSynFile] = fileparts(ctrlConnFile);
ctrlSynToSynFile = sprintf('%s_synToSynDists.mat', ctrlSynToSynFile);
ctrlSynToSynFile = fullfile(fileparts(connFile), ctrlSynToSynFile);
ctrlSynToSyn = load(ctrlSynToSynFile);

%% Prepare data
synT = connectEM.Connectome.buildSynapseTable(conn, syn);
allAxonIds = find(conn.axonMeta.synCount);

plotConfigs = struct;
plotConfigs(1).synIds = find( ...
    synT.isSpine & ismember( ...
    synT.preAggloId, allAxonIds));
plotConfigs(1).title = 'all spine synapses';
plotConfigs(1).tag = 'sp';

plotConfigs(2).synIds = find( ...
    synT.isSpine & ismember( ...
    synT.preAggloId, axonClasses(3).axonIds));
plotConfigs(2).title = 'thalamocortical spine synapses';
plotConfigs(2).tag = 'tc sp';

plotConfigs(3).synIds = find( ...
    synT.isSpine & ismember( ...
    synT.preAggloId, axonClasses(4).axonIds));
plotConfigs(3).title = 'corticocortical spine synapses';
plotConfigs(3).tag = 'cc sp';

ctrlSynT = connectEM.Connectome.buildSynapseTable(ctrlConn, ctrlSyn);
ctrlAllAxonIds = find(ctrlConn.axonMeta.synCount);

ctrlPlotConfigs = struct;
ctrlPlotConfigs(1).synIds = find( ...
    ctrlSynT.isSpine & ismember( ...
    ctrlSynT.preAggloId, ctrlAllAxonIds));
ctrlPlotConfigs(1).title = 'all spine synapses (control)';
ctrlPlotConfigs(1).tag = 'sp (ctrl)';

%% Plot distribution of synapse size
connectEM.Consistency.plotSizeHistogram(info, synT, plotConfigs(1));
connectEM.Consistency.plotSizeHistogram(info, synT, plotConfigs(2:3));

%% Plot histogram of degree of coupling
connectEM.Consistency.plotCouplingHistogram( ...
    info, synT, plotConfigs(1), 'normalization', 'count');
connectEM.Consistency.plotCouplingHistogram( ...
    info, ctrlSynT, ctrlPlotConfigs(1), 'normalization', 'count');

%% Synapse areas vs. degree of coupling
clear cur*;
curPlotCouplings = 1:5;
curConfigs = [ ...
    struct('synT', synT, 'plotConfig', plotConfigs(1)), ...
    struct('synT', ctrlSynT, 'plotConfig', ctrlPlotConfigs)];

for curConfig = reshape(curConfigs, 1, [])
    curPlotConfig = curConfig.plotConfig;

    curSynT = curConfig.synT(curPlotConfig.synIds, :);
    [~, ~, curSynT.pairId] = unique(curSynT(:, ...
        {'preAggloId', 'postAggloId'}), 'rows');

    curCouplings = accumarray(curSynT.pairId, 1);
    curSynT.coupling = curCouplings(curSynT.pairId);

    curPlotConfigs = arrayfun( ...
        @(c) struct( ...
            'coupling', c, ...
            'synIds', curPlotConfig.synIds(curSynT.coupling == c), ...
            'title', sprintf('%d-fold %s', c, curPlotConfig.title)), ...
        curPlotCouplings);
    
    connectEM.Consistency.plotSizeHistogram( ...
        info, curConfig.synT, curPlotConfigs, ...
        'scale', 'log', 'title', curPlotConfig.title);

    fig = ...
        connectEM.Consistency.plotSizeBoxPlot( ...
            info, curConfig.synT, curPlotConfigs, ...
            'title', curPlotConfig.title);
    fig.Position(3:4) = [250, 420];
end

%% Illustrate synapse size similarity
clear cur*;
curPlotConfig = plotConfigs(1);
curPairConfigs = ...
    connectEM.Consistency.buildPairConfigs(synT, curPlotConfig);

for curPairConfig = curPairConfigs(1:(end - 1))
    curPairConfig(2) = curPairConfigs(end); %#ok
    curFig = connectEM.Consistency.plotVariabilityPaired( ...
        info, synT, curPlotConfig, curPairConfig, 'lineCount', 10);
    curFig.Position(3:4) = [700, 330];
end

%% Synapse area variability
clear cur*;
curConfigs = [ ...
    struct('synT', synT, 'plotConfigs', plotConfigs), ...
    struct('synT', ctrlSynT, 'plotConfigs', ctrlPlotConfigs)];

for curConfig = curConfigs
    for curPlotConfig = curConfig.plotConfigs
        curPairConfigs = ...
            connectEM.Consistency.buildPairConfigs( ...
                curConfig.synT, curPlotConfig);

        curFig = ...
            connectEM.Consistency.plotVariabilityHistogram( ...
                info, curConfig.synT, curPlotConfig, curPairConfigs(:));
        curFig.Position(3:4) = [370, 540];

       [curLearnedFrac, curUnlearnedFrac, curCvThresh] = ...
            connectEM.Consistency.calculateLearnedFraction( ...
                curConfig.synT, curPairConfigs(1), curPairConfigs(end));

        fprintf('%s\n', curPlotConfig.title);
        fprintf('→ Learned fraction: %.1f %%\n', 100 * curLearnedFrac);
        fprintf('→ Unlearned fraction: %.1f %%\n', 100 * curUnlearnedFrac);
        fprintf('→ CV threshold: %.2f\n', curCvThresh);
        fprintf('\n');

        fprintf('* Significance tests\n');
        fprintf('  p-values for unexpected synapse size similarity\n');
        for curPairConfig = curPairConfigs(1:(end - 1))
            curPValue = ...
                connectEM.Consistency.testVariability( ...
                    curConfig.synT, curPairConfig, curPairConfigs(end));
            fprintf('→ %s: %g\n', curPairConfig.title, curPValue);
        end

        fprintf('\n');
    end
end

%% Synapse size variability vs. degrees of coupling
clear cur*;
curPlotCouplings = 2:5;
curConfigs = [ ...
    struct('synT', synT, 'plotConfig', plotConfigs(1)), ...
    struct('synT', ctrlSynT, 'plotConfig', ctrlPlotConfigs)];

for curConfig = curConfigs
    curPlotConfig = curConfig.plotConfig;

    curSynT = curConfig.synT(curPlotConfig.synIds, :);
    curSynT.synId = curPlotConfig.synIds;
    
   [~, ~, curSynT.pairId] = unique(curSynT(:, ...
        {'preAggloId', 'postAggloId'}), 'rows');
    curSynT = sortrows(curSynT, 'pairId');

    curCouplings = accumarray(curSynT.pairId, 1);
    curSynT.coupling = curCouplings(curSynT.pairId);
    
    curPairs = arrayfun( ...
        @(c) struct( ...
            'synIdPairs', transpose(reshape( ...
                curSynT.synId(curSynT.coupling == c), c, [])), ...
            'title', sprintf( ...
                '%d-fold %s', c, curPlotConfig.title)), ...
        curPlotCouplings);
    curCtrlPairs = arrayfun( ...
        @(c) struct( ...
            'synIdPairs', reshape(curPlotConfig.synIds(randperm( ...
                c * floor(numel(curPlotConfig.synIds) / c))), [], c), ...
            'title', sprintf( ...
                '%d-fold %s (random)', c, curPlotConfig.title)), ...
        curPlotCouplings);
    
    curPlotPairs = cat(1, curPairs, curCtrlPairs);
    curPlotConfigs = repelem(curPlotConfig, 1, numel(curPlotCouplings));
    
    curFig = ...
        connectEM.Consistency.plotVariabilityHistogram( ...
            info, curConfig.synT, curPlotConfigs, curPlotPairs);
    curFig.Position(3:4) = [1300, 500];
    
    fprintf('* Significance tests\n');
    fprintf('  p-values for unexpected synapse size similarity\n');
    for curPairConfig = curPlotPairs
        curPValue = ...
            connectEM.Consistency.testVariability( ...
                curConfig.synT, curPairConfig(1), curPairConfig(2));
        fprintf('→ %s: %g\n', curPairConfig(1).title, curPValue);
    end
    
    fprintf('\n');
end
            
%% Variability of largest two synapses
curPlotCouplings = 2:5;
[~, curCouplings, curPlotConfigs] = ...
    ndgrid(1, curPlotCouplings, plotConfigs);

curPairConfigs = @(conf, coup) ...
    connectEM.Consistency.buildLargestPairConfigs(synT, conf, coup);
curPairConfigs = cell2mat(arrayfun( ...
    @(conf, coup) reshape(curPairConfigs(conf, coup), [], 1), ...
    curPlotConfigs, curCouplings, 'UniformOutput', false));

curTitles = arrayfun( ...
    @(coup, conf) sprintf('%d-fold %s', coup, conf.title), ...
    curCouplings, curPlotConfigs, 'UniformOutput', false);
[curPlotConfigs.title] = deal(curTitles{:});

connectEM.Consistency.plotVariabilityHistogram( ...
    info, synT, curPlotConfigs, curPairConfigs);

%% Variability vs. distance
clear cur*;
curConfigs = [ ...
    struct( ...
        'synT', synT, 'plotConfig', plotConfigs(1), ...
        'synToSyn', synToSyn, 'maxDistUm', []), ...
    struct( ...
        'synT', synT, 'plotConfig', plotConfigs(1), ...
        'synToSyn', synToSyn, 'maxDistUm', 20), ...
    struct( ...
        'synT', ctrlSynT, 'plotConfig', ctrlPlotConfigs, ...
        'synToSyn', ctrlSynToSyn, 'maxDistUm', []), ...
    struct( ...
        'synT', ctrlSynT, 'plotConfig', ctrlPlotConfigs, ...
        'synToSyn', ctrlSynToSyn, 'maxDistUm', 20)];

for curConfig = curConfigs
    curSynT = curConfig.synT;
    curPlotConfig = curConfig.plotConfig;
    curSynToSyn = curConfig.synToSyn;
    
    curPairConfig = ...
        connectEM.Consistency.buildPairConfigs(curSynT, curPlotConfig);
    curPairConfig = curPairConfig(1);
    
    connectEM.Consistency.plotVariabilityVsDistance( ...
        curSynT, curSynToSyn, curPairConfig, ...
        'maxDistUm', curConfig.maxDistUm);
end
