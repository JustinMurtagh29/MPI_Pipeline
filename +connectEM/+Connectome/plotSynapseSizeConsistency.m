% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, syn, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);

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

%% Plot distribution of synapse size
connectEM.Consistency.plotSizeHistogram(info, synT, plotConfigs(1));
connectEM.Consistency.plotSizeHistogram(info, synT, plotConfigs(2:3));

%% Plot histogram of degree of coupling
connectEM.Consistency.plotCouplingHistogram( ...
    info, synT, plotConfigs(1), 'normalization', 'count');
connectEM.Consistency.plotCouplingHistogram( ...
    info, synT, plotConfigs(2:3), 'normalization', 'probability');

%% Synapse areas vs. degree of coupling
curPlotCouplings = 1:5;
curPlotConfig = plotConfigs(1);

curSynT = synT(curPlotConfig.synIds, :);
[~, ~, curSynT.pairId] = unique(curSynT(:, ...
    {'preAggloId', 'postAggloId'}), 'rows');

curCouplings = accumarray(curSynT.pairId, 1);
curSynT.coupling = curCouplings(curSynT.pairId);

curPlotConfigs = arrayfun( ...
    @(c) struct( ...
        'synIds', curPlotConfig.synIds(curSynT.coupling == c), ...
        'title', sprintf('%d-fold %s', c, curPlotConfig.title)), ...
	curPlotCouplings);

connectEM.Consistency.plotSizeHistogram(info, synT, curPlotConfigs);

%% Synapse area variability
for curPlotConfig = plotConfigs
    curPairConfigs = ...
        connectEM.Consistency.buildPairConfigs(synT, curPlotConfig);
    connectEM.Consistency.plotVariabilityHistogram( ...
        info, synT, curPlotConfig, curPairConfigs);
end

%% Variability of largest two synapses
curPlotCouplings = 2:5;

for curPlotConfig = plotConfigs
    curPlotConfig.rawTitle = curPlotConfig.title;
    
    for curCoupling = curPlotCouplings
        curPairConfigs = ...
            connectEM.Consistency.buildLargestPairConfigs( ...
                synT, curPlotConfig, curCoupling);
        curPlotConfig.title = sprintf( ...
            '%d-fold %s', curCoupling, curPlotConfig.rawTitle);
        connectEM.Consistency.plotVariabilityHistogram( ...
            info, synT, curPlotConfig, curPairConfigs);
    end
end

%% ASI sizes (in log scale)
plotCouplings = 2:5;
fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [960, 1090];

for curCouplingIdx = 1:numel(plotCouplings)
    curCoupling = plotCouplings(curCouplingIdx);
    
    synAreas = zeros(0, 1);
    synGroups = zeros(0, 1);

    for curClassIdx = 1:numel(plotConfigs)
        curAxonIds = plotConfigs(curClassIdx).axonIds;
        curSynIds = plotConfigs(curClassIdx).synIds;
        
        curSynGroups = synT(curSynIds, {'preAggloId', 'postAggloId'});
       [~, ~, curSynGroups] = unique(curSynGroups, 'rows');
       
        curSynAreas = accumarray( ...
            curSynGroups, synT.area(curSynIds), ...
            [], @(areas) {sort(areas, 'descend')});
        
        curSynGroups = accumarray(curSynGroups, 1);
        curSynAreas(curSynGroups ~= curCoupling) = [];
        
        curSynGroups = transpose(repmat( ...
            1:curCoupling, 1, numel(curSynAreas)));
        curSynGroups = (curSynGroups - 1) ...
            * numel(plotConfigs) + curClassIdx;
        curSynAreas = cell2mat(curSynAreas);
        
        synAreas = [synAreas; curSynAreas(:)]; %#ok
        synGroups = [synGroups; curSynGroups(:)]; %#ok
    end
    
    curGroupLabels = arrayfun(@(i) sprintf( ...
        '(%d) (%s)', ...
        ceil(i / numel(plotConfigs)), ...
        plotConfigs(1 + mod(i - 1, numel(plotConfigs))).tag), ...
        unique(synGroups), 'UniformOutput', false);
    curGroupIds = mod(synGroups - 1, numel(plotConfigs));
    
    figure(fig);
    curAx = subplot(numel(plotCouplings), 1, curCouplingIdx);
    synAreasLog = log10(synAreas);
    
    boxplot( ...
        curAx, synAreasLog, synGroups, ...
        'Labels', curGroupLabels, ...
        'Colors', curAx.ColorOrder, ...
        'ColorGroup', curGroupIds, ...
        'Symbol', '.');
    
    curAx.TickDir = 'out';
    curAx.XTickLabelRotation = 30;
end

xlabel(curAx, 'Synapse');
ylabel(curAx, 'log_{10}(ASI)');

xMax = max(arrayfun(@(a) a.XLim(end), fig.Children));
yLims = cell2mat(arrayfun( ...
    @(a) a.YLim, fig.Children(:), ...
    'UniformOutput', false));
yLims = prctile(yLims(:), [0, 100]);

[fig.Children.XLim] = deal([0, xMax]);
[fig.Children.YLim] = deal(yLims);

annotation( ...
    fig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
	'String', {info.filename; info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');

%% ASI variability for all combinations
plotCouplings = 2:5;
cvOf = @(d) std(d, 0, 2) ./ mean(d, 2);

fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [960, 1090];

for curCouplingIdx = 1:numel(plotCouplings)
    curCoupling = plotCouplings(curCouplingIdx);
    
    curPairs = sortrows(combnk(1:curCoupling, 2));
   [~, curSortIds] = sort(diff(curPairs, 1, 2), 'ascend');
    curPairs = curPairs(curSortIds, :);
    
    cvVals = zeros(0, 1);
    cvGroups = zeros(0, 1);

    for curClassIdx = 1:numel(plotConfigs)
        curAxonIds = plotConfigs(curClassIdx).axonIds;
        curSynIds = plotConfigs(curClassIdx).synIds;

       [~, ~, curCvGroups] = unique(synT( ...
            curSynIds, {'preAggloId', 'postAggloId'}), 'rows');
        curSynAreas = accumarray( ...
            curCvGroups, synT.area(curSynIds), ...
            [], @(areas) {sort(areas, 'descend')});
        
        curCvGroups = accumarray(curCvGroups, 1);
        curSynAreas(curCvGroups ~= curCoupling) = [];
        
        curCvVals = cellfun( ...
            @(asis) cvOf(reshape(asis(curPairs), [], 2))', ...
            curSynAreas, 'UniformOutput', false);
        curCvVals = cell2mat(curCvVals);
        
        curCvGroups = repelem( ...
            1:size(curPairs, 1), size(curCvVals, 1));
        curCvGroups = (curCvGroups - 1) ...
            * numel(plotConfigs) + curClassIdx;
        
        cvVals = [cvVals; curCvVals(:)]; %#ok
        cvGroups = [cvGroups; curCvGroups(:)]; %#ok
    end
    
    curGroupLabels = arrayfun(@(i) sprintf( ...
        '(%d, %d) (%s)', ...
        curPairs(ceil(i / numel(plotConfigs)), :), ...
        plotConfigs(1 + mod(i - 1, numel(plotConfigs))).tag), ...
        unique(cvGroups), 'UniformOutput', false);
    curGroupIds = mod(cvGroups - 1, numel(plotConfigs));
    
    figure(fig);
    curAx = subplot(numel(plotCouplings), 1, curCouplingIdx);
    
    boxplot( ...
        curAx, cvVals, cvGroups, ...
        'Labels', curGroupLabels, ...
        'Colors', curAx.ColorOrder, ...
        'ColorGroup', curGroupIds, ...
        'Symbol', '.');
    
    curAx.TickDir = 'out';
    curAx.XTickLabelRotation = 30;
end

xlabel(curAx, 'Synapse pair');
ylabel(curAx, 'Synapse variability (CV)');

xMax = max(arrayfun(@(a) a.XLim(end), fig.Children));
yMax = max(arrayfun(@(a) a.YLim(end), fig.Children));

[fig.Children.XLim] = deal([0, xMax]);
[fig.Children.YLim] = deal([0, yMax]);

annotation( ...
    fig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
	'String', { ...
        'Variability for all synapse combinations'; ...
        info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');
