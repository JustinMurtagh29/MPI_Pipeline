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
curPlotConfig = plotConfigs(1);

curSynT = synT;
curSynT.id = reshape(1:size(curSynT, 1), [], 1);
curSynT = curSynT(curPlotConfig.synIds, :);

[~, curUniSynT, curSynT.pairId] = unique( ...
    curSynT(:, {'preAggloId', 'postAggloId'}), 'rows');
curUniSynT = curSynT(curUniSynT, :);

% Same-axon same-dendrite
curSaSdIds = curSynT;
curSaSdIds.coupling = curCouplings(curSaSdIds.pairId);
curSaSdIds(curSaSdIds.coupling ~= 2, :) = [];
curSaSdIds = reshape(curSaSdIds.id, 2, [])';

% Same-axon different-dendrite
curSaDdIds = cell2mat(accumarray( ...
    curUniSynT.preAggloId, curUniSynT.id, [], @(ids) { ...
    reshape(ids(1:(2 * floor(numel(ids) / 2))), [], 2)}));

% Different-axon same-dendrite
curDaSdIds = cell2mat(accumarray( ...
    curUniSynT.postAggloId, curUniSynT.id, [], @(ids) { ...
    reshape(ids(1:(2 * floor(numel(ids) / 2))), [], 2)}));

% Random pairs
rng(0);
curRandIds = randperm(2 * floor(size(curSynT, 1) / 2));
curRandIds = reshape(curSynT.id(curRandIds), [], 2);

curPlotConfigs = struct;
curPlotConfigs(1).synPairIds = curSaSdIds;
curPlotConfigs(1).title = 'Same-axon same-dendrite';

curPlotConfigs(2).synPairIds = curSaDdIds;
curPlotConfigs(2).title = 'Same-axon different-dendrite';

curPlotConfigs(3).synPairIds = curDaSdIds;
curPlotConfigs(3).title = 'Same-axon different-dendrite';

curPlotConfigs(4).synPairIds = curRandIds;
curPlotConfigs(4).title = 'Random pairs';

%% Synapse area variability
cvVals = zeros(0, 1);
cvGroups = zeros(0, 1);

for curClassIdx = 1:numel(plotConfigs)
    curAxonIds = plotConfigs(curClassIdx).axonIds;
    curSynIds = plotConfigs(curClassIdx).synIds;
    
   [~, ~, curCvGroups] = unique(synT( ...
        curSynIds, {'preAggloId', 'postAggloId'}), 'rows');
    curCvVals = accumarray( ...
        curCvGroups, synT.area(curSynIds), [], ...
        @(areas) std(areas) / mean(areas));
    curCvGroups = accumarray(curCvGroups, 1);

    % remove single-spine case
    curCvVals(curCvGroups < 2) = [];
    curCvGroups(curCvGroups < 2) = [];
    
    % translate group id
    curCvGroups = (curCvGroups - 2) ...
        * numel(plotConfigs) + curClassIdx;
    
    cvVals = [cvVals; curCvVals]; %#ok
    cvGroups = [cvGroups; curCvGroups]; %#ok
end

groupLabels = arrayfun(@(i) sprintf( ...
    '%d (%s)', 1 + ceil(i / numel(plotConfigs)), ...
    plotConfigs(1 + mod(i - 1, numel(plotConfigs))).tag), ...
    unique(cvGroups), 'UniformOutput', false);
groupIds = mod(cvGroups - 1, numel(plotConfigs));

% plot
fig = figure();
fig.Color = 'white';
ax = axes(fig);

boxplot( ...
    ax, cvVals, cvGroups, ...
    'Labels', groupLabels, ...
    'Colors', ax.ColorOrder, ...
    'ColorGroup', groupIds, ...
    'Symbol', '.');

ylabel(ax, 'Variability of synapse areas (CV)');
xlabel(ax, 'Synapses per connection');

title( ...
    ax, {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);

ax.YLim(1) = 0;
ax.TickDir = 'out';

%% do comparisons against null hypothesis
controlCouplings = 2:5;
cvOf = @(d) std(d, 0, 2) ./ mean(d, 2);

for curCoupling = controlCouplings
    curFig = figure();
    curFig.Color = 'white';
    curFig.Position(3:4) = [520, 530];
    
    for curClassIdx = 1:numel(plotConfigs)
        curAxonIds = plotConfigs(curClassIdx).axonIds;
        curSynIds = plotConfigs(curClassIdx).synIds;
        
        % actual CVs
       [~, ~, curSynCoupling] = unique(synT( ...
            curSynIds, {'preAggloId', 'postAggloId'}), 'rows');
        curSynAreas = accumarray( ...
            curSynCoupling, synT.area(curSynIds), [], ...
            @(a) {reshape(sort(a, 'descend'), 1, [])});
        
        curSynCoupling = accumarray(curSynCoupling, 1);
        curSynAreas(curSynCoupling ~= curCoupling) = [];
        
        curSynAreas = cell2mat(curSynAreas);
        curCvs = cvOf(curSynAreas(:, 1:2));
        
        % control
        rng(0);
        curCtrlCount = curCoupling * ...
            floor(numel(curSynIds) / curCoupling);

        curCtrlVals = randperm(numel(curSynIds), curCtrlCount);
        curCtrlVals = reshape(curSynIds(curCtrlVals), [], curCoupling);
        curCtrlVals = synT.area(curCtrlVals);

        curCtrlVals = sort(curCtrlVals, 2, 'descend');
        curCtrlVals = cvOf(curCtrlVals(:, 1:2));
        
        curAx = subplot(numel(plotConfigs), 1, curClassIdx);
        axis(curAx, 'square');
        hold(curAx, 'on');
        
        histogram( ...
            curAx, curCtrlVals, linspace(0, 1.5, 21), ...
            'Normalization', 'probability', ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2, ...
            'FaceAlpha', 1);
        histogram( ...
            curAx, curCvs, linspace(0, 1.5, 21), ...
            'Normalization', 'probability', ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2, ...
            'FaceAlpha', 1);
        
        title( ...
            curAx, plotConfigs(curClassIdx).title, ...
            'FontWeight', 'normal', 'FontSize', 10);
        curAx.TickDir = 'out';
    end
    
    annotation( ...
        curFig, ...
        'textbox', [0, 0.9, 1, 0.1], ...
        'String', { ...
            sprintf('%d-fold coupled neurites', curCoupling);
            info.git_repos{1}.hash}, ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    
    legend(curAx, 'Control', 'Observed');
    xlabel(curAx, 'CV of largest two AIS');
    ylabel(curAx, 'Probability');
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
