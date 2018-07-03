% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-linearized_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto.mat');

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, syn, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);

[~, synToSynFile] = fileparts(connFile);
synToSynFile = sprintf('%s_synToSynDists.mat', synToSynFile);
synToSynFile = fullfile(fileparts(connFile), synToSynFile);
synToSyn = load(synToSynFile);

% Loading spine head agglomerates
shAgglos = load(shFile, 'shAgglos');
shAgglos = shAgglos.shAgglos;

% Loading augmented graph
graph = Graph.load(rootDir);
graph(~graph.borderIdx, :) = [];

borderAreas = fullfile(rootDir, 'globalBorder.mat');
borderAreas = load(borderAreas, 'borderArea2');
borderAreas = borderAreas.borderArea2;

graph.borderArea = borderAreas(graph.borderIdx);
clear borderAreas;

% HACK(amotta): For some reason there exist borders, for which
% `physicalBorderArea2` is zero. This seems wrong.
%   In order not to be affected by this issue, let's set the area of these
% borders to NaN. This will result in a total axon-spine interface area of
% NaN, which we can remove by brute force later on.
%
% Corresponding issue on GitLab:
% https://gitlab.mpcdf.mpg.de/connectomics/auxiliaryMethods/issues/16
graph.borderArea(~graph.borderArea) = nan;
graph(:, {'prob', 'borderIdx'}) = [];

%% Build axon-spine interface areas
asiT = ...
    connectEM.Connectome.buildAxonSpineInterfaces( ...
        param, graph, shAgglos, conn, syn);
asiT(asiT.type ~= 'PrimarySpine', :) = [];

% HACK(amotta): This is the counter-part to the above hack. (This is not
% really needed, since `buildAxonSpineInterface` does exactly the same
% thing internally. But let's be explicit and future-proof.)
asiT(isnan(asiT.area), :) = [];

synT = asiT;
synT.isSpine(:) = true;

%% Prepare data
allAxonIds = unique(asiT.preAggloId);

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

%% Report synapse sizes
clear cur*;
fprintf('Synapse sizes\n');
for curConfig = plotConfigs
    curSynT = synT(curConfig.synIds, :);
    curLogAsi = log10(curSynT.area);
    curMeanLog10Asi = mean(curLogAsi);
    curStdLog10Asi = std(curLogAsi);
    
    fprintf( ...
        '* log10(%s ASI): %f ± %f (mean ± std)\n', ...
        curConfig.title, curMeanLog10Asi, curStdLog10Asi)
end

%% Plot distribution of synapse size
connectEM.Consistency.plotSizeHistogram( ...
    info, synT, plotConfigs(1), 'scale', 'log');
connectEM.Consistency.plotSizeHistogram( ...
    info, synT, plotConfigs(2:3), 'scale', 'log');

%% Report number of occurences for degree of coupling
clear cur*;
curPlotConfig = plotConfigs(1);
curSynT = synT(curPlotConfig.synIds, :);
[~, ~, curDegreeOfCoupling] = unique(curSynT( ...
    :, {'preAggloId', 'postAggloId'}), 'rows');
curDegreeOfCoupling = accumarray(curDegreeOfCoupling, 1);

docT = table;
[docT.degreeOfCoupling, ~, curDocCount] = unique(curDegreeOfCoupling);

docT.occurences = accumarray(curDocCount, 1);
docT(docT.degreeOfCoupling == 1, :) = [];

disp(docT)

%% Plot histogram of degree of coupling
connectEM.Consistency.plotCouplingHistogram( ...
    info, synT, plotConfigs(1), 'normalization', 'count');

%% Synapse areas vs. degree of coupling
clear cur*;
curPlotCouplings = 1:4;
curConfigs = struct('synT', synT, 'plotConfig', plotConfigs(1));

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

   [curFig, curFit] = ...
        connectEM.Consistency.plotSizeBoxPlot( ...
            info, curConfig.synT, curPlotConfigs, ...
            'title', curPlotConfig.title);
    curFig.Position(3:4) = [250, 420];
    
    curAsiGainVsDoc = 10 ^ curFit.p1;
    fprintf('Fit results for "%s"\n', curPlotConfig.title);
    fprintf('* Fold ASI per degree of coupling: %f\n', curAsiGainVsDoc);
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
curConfigs = struct('synT', synT, 'plotConfigs', plotConfigs);

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
        curOpenFrac = 1 - curLearnedFrac - curUnlearnedFrac;
        
        % Quantitative reporting
        fprintf('**%s**\n', curPlotConfig.title);
        curSynIdPairs = curPairConfigs(1).synIdPairs;
        
        curSynT = curConfig.synT( ...
            curSynIdPairs(:, 1), {'preAggloId', 'postAggloId'});
        curSynT.areas = curConfig.synT.area(curSynIdPairs);
        curSynT.cv = std(curSynT.areas, 0, 2) ./ mean(curSynT.areas, 2);
        curSynT.low = curSynT.cv < curCvThresh;
        
        curAxonT = curSynT(:, {'preAggloId'});
       [curAxonT, ~, curGroupIds] = unique(curAxonT);
        curAxonT.allConn = accumarray(curGroupIds, 1);
        curAxonT.lowConn = accumarray(curGroupIds, curSynT.low);
        curAxonT = sortrows(curAxonT, 'lowConn', 'descend');
        
        curLowThresh = round(height(curSynT) * curLearnedFrac);
        curAxonFracLow = find(cumsum(curAxonT.lowConn) > curLowThresh, 1);
        curAxonFracLow = curAxonFracLow / height(curAxonT) %#ok
        curAxonFracUp = curLowThresh / height(curAxonT) %#ok
        
        fprintf('CV between pairs (mean ± std)\n');
        for curPairConfig = curPairConfigs
            curCvs = curConfig.synT.area(curPairConfig.synIdPairs);
            curCvs = std(curCvs, 0, 2) ./ mean(curCvs, 2);
            
            fprintf( ...
                '* %s: %f ± %f\n', ...
                curPairConfig.title, ...
                mean(curCvs), std(curCvs));
        end
            
        fprintf('\n');
        
        fprintf('* Learned fraction: %.1f %%\n', 100 * curLearnedFrac);
        fprintf('* Possibly learned fraction: %.1f %%\n', 100 * curOpenFrac);
        fprintf('* Unlearned fraction: %.1f %%\n', 100 * curUnlearnedFrac);
        fprintf('* CV threshold: %.2f\n', curCvThresh);
        fprintf('\n');

        fprintf('Significance tests\n');
        fprintf('p-values for unexpected synapse size similarity\n');
        for curPairConfig = curPairConfigs(1:(end - 1))
            curPValue = ...
                connectEM.Consistency.testVariability( ...
                    curConfig.synT, curPairConfig, curPairConfigs(end));
            fprintf('* %s: %g\n', curPairConfig.title, curPValue);
        end
        
        fprintf('\n');
    end
end

%% Correlation of synapse size correlation with synapse size
clear cur*;
curMinCv = 0.0;
curMaxCv = 0.25;

curSynT = synT;
curConfig = plotConfigs(1);

curPairConfigs = ...
    connectEM.Consistency.buildPairConfigs(curSynT, curConfig);

curCtrlPairConfig = ...
    connectEM.Consistency.buildPairConfigs(curSynT, struct( ...
        'synIds', curPairConfigs(1).synIdPairs(:)));
curCtrlPairConfig = curCtrlPairConfig(end);
curCtrlPairConfig.title = sprintf( ...
    'Random pairs from above set (n = %d)', ...
    size(curCtrlPairConfig.synIdPairs, 1));

curFig = ...
	connectEM.Consistency.plotVariabilityHistogram( ...
        info, curSynT, curConfig, [curPairConfigs(1); curCtrlPairConfig]);
curFig.Position(3:4) = [370, 540];

curPairT = table;
curPairT.areas = curSynT.area(curPairConfigs(1).synIdPairs);
curPairT.areas = sort(curPairT.areas, 2);
curPairT.cv = std(curPairT.areas, 0, 2) ./ mean(curPairT.areas, 2);

curRandT = table;
curRandT.areas = curSynT.area(curPairConfigs(end).synIdPairs);
curRandT.areas = sort(curRandT.areas, 2);
curRandT.cv = std(curRandT.areas, 0, 2) ./ mean(curRandT.areas, 2);

curCtrlT = table;
curCtrlT.areas = curSynT.area(curCtrlPairConfig.synIdPairs);
curCtrlT.areas = sort(curCtrlT.areas, 2);
curCtrlT.cv = std(curCtrlT.areas, 0, 2) ./ mean(curCtrlT.areas, 2);

curFig = figure;
curFig.Color = 'white';
curFig.Position(3:4) = [420, 820];

curAx = subplot(2, 1, 1);
hold(curAx, 'on');

curBinEdges = linspace(-1.5, 0.5, 11);
curHist = @(t) histogram(curAx, ...
    mean(log10(t.areas( ...
        t.cv > curMinCv ...
      & t.cv < curMaxCv, :)), 2), ...
	'BinEdges', curBinEdges, ...
    'DisplayStyle', 'stairs', ...
    'LineWidth', 2);
curPairHist = curHist(curPairT);
curCtrlHist = curHist(curCtrlT);

xlabel(curAx, 'log_{10}(ASI area [µm²])');
ylabel(curAx, 'Probability');

curAx = subplot(2, 1, 2);
histogram(curAx, ...
    'BinCounts', max( ...
        curPairHist.BinCounts ...
      - curCtrlHist.BinCounts, 0), ...
	'BinEdges', curBinEdges, ...
    'DisplayStyle', 'stairs', ...
    'LineWidth', 2);

xlabel(curAx, 'log_{10}(ASI area [µm²])');
ylabel(curAx, 'Observed vs. expected synapse pairs');

curAxes = flip(cat(1, curFig.Children));
set(curAxes, ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'XLim', curBinEdges([1, end]), ...
    'PlotBoxAspectRatio', [1, 1, 1],...
	'DataAspectRatioMode', 'auto');

curAx = curAxes(1);
curPos = curAx.Position;
curLeg = legend(curAx, { ...
    'Same-axon same-dendrite spine synapses', ...
    'Random pairs from above set of synapses'}, ...
    'Location', 'SouthOutside');
curLeg.Box = 'off';
curAx.Position = curPos;

annotation( ...
    curFig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
	'String', { ...
        info.filename; info.git_repos{1}.hash; ....
        sprintf('Pairs with %g < CV < %g', curMinCv, curMaxCv)}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');

%% Synapse size variability vs. degrees of coupling
clear cur*;
curPlotCouplings = 2:5;
curConfigs = struct('synT', synT, 'plotConfig', plotConfigs(1));

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
        'synToSyn', synToSyn, 'maxDistUm', 20)];
        
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
