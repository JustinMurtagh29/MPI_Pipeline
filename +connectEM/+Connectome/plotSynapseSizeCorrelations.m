% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
synFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v3_ax_spine_clustered_classified.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_b_linearized_ax_spine_syn_clust.mat');

info = Util.runInfo();

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[interSynDir, interSynFile] = fileparts(connFile);
interSynFile = sprintf('%s_intersynapse.mat', interSynFile);
interSynFile = fullfile(interSynDir, interSynFile);
clear interSynDir;

interSyn = load(interSynFile);

conn = load(connFile);
syn = load(synFile);

%% limit synapses
synT = connectEM.Connectome.buildSynapseTable(conn, syn);
synT.isSoma = (conn.denMeta.targetClass(synT.postAggloId) == 'Somata');
synAreaLim = prctile(synT.area, 99.9);

%% define axon classes
axonClasses = struct;

% spine vs. shaft synapses
axonClasses(1).axonIds = find(conn.axonMeta.synCount);
axonClasses(1).synIds = find(synT.isSpine);
axonClasses(1).title = 'Spine synapses';
axonClasses(1).tag = 'Sp';

%% plot distribution of synapse size
binEdges = linspace(0, synAreaLim, 51);

fig = figure();
fig.Color = 'white';
ax = axes(fig);

hold(ax, 'on');
ax.TickDir = 'out';

for curClassIdx = 1:numel(axonClasses)
    curAxonIds = axonClasses(curClassIdx).axonIds;
    curSynIds = axonClasses(curClassIdx).synIds;
    
    histogram( ...
        ax, synT.area(curSynIds), ...
        'BinEdges', binEdges, ...
        'Normalization', 'probability', ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2);
end

xlim(ax, binEdges([1, end]));
xlabel(ax, 'Synapse area (µm²)');
ylabel(ax, 'Probability');

title( ...
   {'Synapse area distribution'; info.git_repos{1}.hash}, ...
    'FontWeight', 'norma', 'FontSize', 10);
legend(ax, {axonClasses.title}, 'Location', 'NorthEast');
    
fig.Position(3:4) = [820, 475];

%% plot histogram over no. of synapse per neurite pair
fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [495, 400];

ax = axes(fig);
axis(ax, 'square');
hold(ax, 'on');

ax.YScale = 'log';
ax.TickDir = 'out';

for curClassIdx = 1:numel(axonClasses)
    curAxonIds = axonClasses(curClassIdx).axonIds;
    curSynIds = axonClasses(curClassIdx).synIds;
    
   [~, ~, curCoupling] = unique(synT( ...
        curSynIds, {'preAggloId', 'postAggloId'}), 'rows');
    curCoupling = accumarray(curCoupling, 1);

    histogram( ...
        ax, curCoupling, ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2, ...
        'FaceAlpha', 1);
end

ylabel(ax, 'Probability');
yticklabels(ax, arrayfun(@num2str, ...
    yticks(ax), 'UniformOutput', false));

xMax = max(arrayfun( ...
@(h) h.BinEdges(end), ax.Children));
xlim(ax, [0.5, xMax]);
xticks(ax, 1:(xMax - 0.5));
xlabel(ax, 'Synapses per connection');

title( ...
    ax, {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
legend( ...
    ax, {axonClasses.title}, ...
    'Location', 'NorthEast', 'Box', 'off');

%% Synapse areas vs. degree of coupling
asiAreas = zeros(0, 1);
asiGroups = zeros(0, 1);

for curClassIdx = 1:numel(axonClasses)
    curAxonIds = axonClasses(curClassIdx).axonIds;
    curSynIds = axonClasses(curClassIdx).synIds;
    
   [~, ~, curAsiGroups] = unique(synT( ...
        curSynIds, {'preAggloId', 'postAggloId'}), 'rows');
    curAsiAreas = accumarray( ...
        curAsiGroups, synT.area(curSynIds), [], @(areas) {areas});
    curAsiGroups = accumarray(curAsiGroups, 1);
    
    curAsiGroups = repelem( ...
        curAsiGroups, cellfun(@numel, curAsiAreas));
    curAsiAreas = cell2mat(curAsiAreas);
    
    % translate group id
    curAsiGroups = (curAsiGroups - 1) ...
        * numel(axonClasses) + curClassIdx;
    
    asiAreas = [asiAreas; curAsiAreas]; %#ok
    asiGroups = [asiGroups; curAsiGroups]; %#ok
end

groupLabels = arrayfun(@(i) sprintf( ...
    '%d (%s)', ceil(i / numel(axonClasses)), ...
    axonClasses(1 + mod(i - 1, numel(axonClasses))).tag), ...
    unique(asiGroups), 'UniformOutput', false);
groupIds = mod(asiGroups - 1, numel(axonClasses));

% plot
fig = figure();
fig.Color = 'white';
ax = axes(fig);

boxplot( ...
    ax, asiAreas, asiGroups, ...
    'Labels', groupLabels, ....
    'Colors', ax.ColorOrder, ...
    'ColorGroup', groupIds, ...
    'Symbol', '.');

xlabel(ax, 'Synapses per connection');
ylabel(ax, 'Synapse area (µm²)');
ylim(ax, [0, synAreaLim]);

title( ...
    ax, {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);

ax.TickDir = 'out';

%% Synapse area variability
cvVals = zeros(0, 1);
cvGroups = zeros(0, 1);

for curClassIdx = 1:numel(axonClasses)
    curAxonIds = axonClasses(curClassIdx).axonIds;
    curSynIds = axonClasses(curClassIdx).synIds;
    
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
        * numel(axonClasses) + curClassIdx;
    
    cvVals = [cvVals; curCvVals]; %#ok
    cvGroups = [cvGroups; curCvGroups]; %#ok
end

groupLabels = arrayfun(@(i) sprintf( ...
    '%d (%s)', 1 + ceil(i / numel(axonClasses)), ...
    axonClasses(1 + mod(i - 1, numel(axonClasses))).tag), ...
    unique(cvGroups), 'UniformOutput', false);
groupIds = mod(cvGroups - 1, numel(axonClasses));

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
    
    for curClassIdx = 1:numel(axonClasses)
        curAxonIds = axonClasses(curClassIdx).axonIds;
        curSynIds = axonClasses(curClassIdx).synIds;
        
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
        
        curAx = subplot(numel(axonClasses), 1, curClassIdx);
        axis(curAx, 'square');
        hold(curAx, 'on');
        
        histogram( ...
            curAx, curCtrlVals, linspace(0, 1.5, 21), ...
            'Normalization', 'probability', ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2);
        histogram( ...
            curAx, curCvs, linspace(0, 1.5, 21), ...
            'Normalization', 'probability', ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2);
        
        title( ...
            curAx, axonClasses(curClassIdx).title, ...
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

    for curClassIdx = 1:numel(axonClasses)
        curAxonIds = axonClasses(curClassIdx).axonIds;
        curSynIds = axonClasses(curClassIdx).synIds;
        
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
            * numel(axonClasses) + curClassIdx;
        curSynAreas = cell2mat(curSynAreas);
        
        synAreas = [synAreas; curSynAreas(:)]; %#ok
        synGroups = [synGroups; curSynGroups(:)]; %#ok
    end
    
    curGroupLabels = arrayfun(@(i) sprintf( ...
        '(%d) (%s)', ...
        ceil(i / numel(axonClasses)), ...
        axonClasses(1 + mod(i - 1, numel(axonClasses))).tag), ...
        unique(synGroups), 'UniformOutput', false);
    curGroupIds = mod(synGroups - 1, numel(axonClasses));
    
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

    for curClassIdx = 1:numel(axonClasses)
        curAxonIds = axonClasses(curClassIdx).axonIds;
        curSynIds = axonClasses(curClassIdx).synIds;

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
            * numel(axonClasses) + curClassIdx;
        
        cvVals = [cvVals; curCvVals(:)]; %#ok
        cvGroups = [cvGroups; curCvGroups(:)]; %#ok
    end
    
    curGroupLabels = arrayfun(@(i) sprintf( ...
        '(%d, %d) (%s)', ...
        curPairs(ceil(i / numel(axonClasses)), :), ...
        axonClasses(1 + mod(i - 1, numel(axonClasses))).tag), ...
        unique(cvGroups), 'UniformOutput', false);
    curGroupIds = mod(cvGroups - 1, numel(axonClasses));
    
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

%% intersynapse distance for all combinations
plotCouplings = 2:5;

[~, ~, neuriteCoupling] = unique( ...
    synT(:, {'preAggloId', 'postAggloId'}), 'rows');
neuriteSyns = accumarray( ...
    neuriteCoupling, (1:size(synT, 1))', [], ...
    @(rows) {sortrows(synT(rows, :), 'area', 'descend')});
neuriteCoupling = accumarray(neuriteCoupling, 1);

fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [1000, 900];

for curCouplingIdx = 1:numel(plotCouplings)
    curCoupling = plotCouplings(curCouplingIdx);
    
    curPairs = sortrows(combnk(1:curCoupling, 2));
   [~, curSortIds] = sort(diff(curPairs, 1, 2), 'ascend');
    curPairs = transpose(curPairs(curSortIds, :));
    
    curPairIds = find(neuriteCoupling == curCoupling);
    curIsds = nan(numel(curPairIds), size(curPairs, 2));
    
    for curPairIdx = 1:numel(curPairIds)
        curPairId = curPairIds(curPairIdx);
        curPairSyns = neuriteSyns{curPairId};
        
        curAxonId = curPairSyns.preAggloId(1);
        curAxonIsds = (interSyn.axonIds == curAxonId);
        
       [~, curSynIds] = ismember( ...
            curPairSyns.id, interSyn.synIds{curAxonIsds});
        curSynIds = curSynIds(curPairs);
        
        curAxonIsds = interSyn.synToSynDists{curAxonIsds};
        curIsds(curPairIdx, :) = arrayfun( ...
            @(one, two) curAxonIsds(one, two), ...
            curSynIds(1, :), curSynIds(2, :));
    end
    
    figure(fig);
    curAx = subplot(numel(plotCouplings), 1, curCouplingIdx);
    boxplot(curAx, curIsds / 1E3);
    
    curAx.TickDir = 'out';
    curAx.XTickLabel = arrayfun( ...
        @(a, b) sprintf('(%d, %d)', a, b), ...
        curPairs(1, :), curPairs(2, :), ...
        'UniformOutput', false);
end

xlabel(curAx, 'Synapse pair');
ylabel(curAx, 'Intersynapse distance (µm)');

xMax = max(arrayfun(@(a) a.XLim(end), fig.Children));
yMax = max(arrayfun(@(a) a.YLim(end), fig.Children));

[fig.Children.XLim] = deal([0, xMax]);
[fig.Children.YLim] = deal([0, yMax]);

annotation( ...
    fig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
	'String', { ...
        'Intersynapse distance all synapse combinations'; ...
        info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');
