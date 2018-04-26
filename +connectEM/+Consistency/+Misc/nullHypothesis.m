% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-classified_spine-syn-clust.mat');

plotSynCounts = 2:5;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

conn = load(connFile);
syn = load(conn.info.param.synFile);

%% Limit synapses
synT = connectEM.Connectome.buildSynapseTable(conn, syn);

%% Define axon and synapse classes
axonClasses = struct;

% spine vs. shaft synapses
axonClasses(1).axonIds = find(conn.axonMeta.synCount);
axonClasses(1).synIds = find(synT.isSpine);
axonClasses(1).title = 'Spine synapses';
axonClasses(1).tag = 'Sp';

%% Plot ASI
axonClass = axonClasses(1);

fig = figure();
fig.Color = 'white';
for curSynCountIdx = 1:numel(plotSynCounts)
    curSynCount = plotSynCounts(curSynCountIdx);
        
    % Observed
   [~, ~, curSynCoupling] = unique(synT( ...
        axonClass.synIds, {'preAggloId', 'postAggloId'}), 'rows');
	curSynAreas = accumarray( ...
        curSynCoupling, synT.area(axonClass.synIds), ...
        [], @(a) {reshape(sort(a, 'descend'), 1, [])});
    
    curSynCoupling = accumarray(curSynCoupling, 1);
    curSynAreas(curSynCoupling ~= curSynCount) = [];
    curSynAreas = cell2mat(curSynAreas);
    
    % Control
    rng(0);
    curCtrlCount = numel(axonClass.synIds);
    curCtrlCount = curCtrlCount - mod(curCtrlCount, curSynCount);
    
    curCtrlAreas = axonClass.synIds(randperm(curCtrlCount));
    curCtrlAreas = reshape(synT.area(curCtrlAreas), [], curSynCount);
    curCtrlAreas = sort(curCtrlAreas, 2, 'descend');
    
    curPlotData = [curSynAreas; curCtrlAreas];
    curPlotData = log10(curPlotData);
    
    curPlotGroups = [ ...
        1 + 2 .* (repmat(1:curSynCount, size(curSynAreas, 1), 1) - 1); ...
        2 + 2 .* (repmat(1:curSynCount, size(curCtrlAreas, 1), 1) - 1)];
    curPlotLabels = arrayfun( ...
        @num2str, 1:curSynCount, ...
        'UniformOutput', false);
    curPlotLabels = [ ...
        strcat(curPlotLabels, {' (obs.)'}); ...
        strcat(curPlotLabels, {' (ctrl.)'})];
    curPlotLabels = curPlotLabels(:);
    
    curAx = subplot(numel(plotSynCounts), 1, curSynCountIdx);
    
    boxplot(curAx, ...
        curPlotData(:), curPlotGroups(:), ...
        'Colors', curAx.ColorOrder(1:2, :));
    
    curAx.TickDir = 'out';
    curAx.XTickLabel = curPlotLabels;
end

[yMin, yMax] = arrayfun(@(c) deal(c.YLim(1), c.YLim(2)), fig.Children);
[fig.Children.YLim] = deal([min(yMin), max(yMax)]);
ylabel(fig.Children(1), 'log_{10}(ASI)');

annotation(fig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'String', {info.filename; info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');

%% Plot ratios
axonClass = axonClasses(1);

fig = figure();
fig.Color = 'white';
for curSynCountIdx = 1:numel(plotSynCounts)
    curSynCount = plotSynCounts(curSynCountIdx);
        
    % Observed
   [~, ~, curSynCoupling] = unique(synT( ...
        axonClass.synIds, {'preAggloId', 'postAggloId'}), 'rows');
	curSynAreas = accumarray( ...
        curSynCoupling, synT.area(axonClass.synIds), ...
        [], @(a) {reshape(sort(a, 'descend'), 1, [])});
    
    curSynCoupling = accumarray(curSynCoupling, 1);
    curSynAreas(curSynCoupling ~= curSynCount) = [];
    curSynAreas = cell2mat(curSynAreas);
    
    curSynRatios = ...
        curSynAreas(:, 2:end) ...
     ./ curSynAreas(:, 1:(end - 1));
    
    % Control
    rng(0);
    curCtrlCount = numel(axonClass.synIds);
    curCtrlCount = curCtrlCount - mod(curCtrlCount, curSynCount);
    
    curCtrlAreas = axonClass.synIds(randperm(curCtrlCount));
    curCtrlAreas = reshape(synT.area(curCtrlAreas), [], curSynCount);
    curCtrlAreas = sort(curCtrlAreas, 2, 'descend');
    
    curCtrlRatios = ...
        curCtrlAreas(:, 2:end) ...
     ./ curCtrlAreas(:, 1:(end - 1));
    
    curPlotData = [curSynRatios; curCtrlRatios];
    
    curPlotGroups = [ ...
        1 + 2 .* (repmat(1:(curSynCount - 1), size(curSynRatios, 1), 1) - 1); ...
        2 + 2 .* (repmat(1:(curSynCount - 1), size(curCtrlRatios, 1), 1) - 1)];
    curPlotLabels = arrayfun( ...
        @num2str, 1:(curSynCount - 1), ...
        'UniformOutput', false);
    curPlotLabels = [ ...
        strcat(curPlotLabels, {' (obs.)'}); ...
        strcat(curPlotLabels, {' (ctrl.)'})];
    curPlotLabels = curPlotLabels(:);
    
    curAx = subplot(numel(plotSynCounts), 1, curSynCountIdx);
    
    boxplot(curAx, ...
        curPlotData(:), curPlotGroups(:), ...
        'Colors', curAx.ColorOrder(1:2, :));
    
    curAx.TickDir = 'out';
    curAx.XTickLabel = curPlotLabels;
end

[fig.Children.YLim] = deal([0, 1]);
ylabel(curAx, 'ASI_{n + 1} / ASI_n');
xlabel(curAx, 'n');

annotation(fig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'String', {info.filename; info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');
