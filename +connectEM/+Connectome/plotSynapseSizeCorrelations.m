% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
synFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v3_ax_spine_clustered.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a_ax_spine_syn_clust.mat');

[interSynDir, interSynFile] = fileparts(connFile);
interSynFile = sprintf('%s_intersynapse.mat', interSynFile);
interSynFile = fullfile(interSynDir, interSynFile);
clear interSynDir;

% set this variable to debug
debugDir = '';

info = Util.runInfo();

%% loading data
% NOTE(amotta): Synapses sizes are contained in the `contactArea` field of 
% `conn.connectomeMeta`. Each cell contains the synapses sizes of the
% correponding entries in `conn.connectome`.
syn = load(synFile);
interSyn = load(interSynFile);

[~, connName] = fileparts(connFile);
conn = connectEM.Connectome.load( ...
    struct('saveFolder', rootDir), connName);

%% limit synapses
synT = connectEM.Connectome.buildSynapseTable(conn, syn);
synT.isSoma = (conn.denMeta.targetClass(synT.postAggloId) == 'Somata');
synAreaLim = prctile(synT.area, 99.9);

%% define axon classes
conn.axonMeta.spineFrac = ...
    conn.axonMeta.spineSynCount ...
    ./ conn.axonMeta.synCount;

axonClasses = struct;

%{
% excitatory vs. thalamocortical
axonClasses(1).axonIds = find( ...
    conn.axonMeta.synCount >= 10 ...
  & conn.axonMeta.spineFrac >= 0.5);
axonClasses(1).synIds = find(ismember( ...
    synT.preAggloId, axonClasses(1).axonIds));
axonClasses(1).title = 'Excitatory axons';
axonClasses(1).tag = 'Exc';

axonClasses(2).axonIds = find(conn.axonMeta.isThalamocortical);
axonClasses(2).synIds = find(ismember( ...
    synT.preAggloId, axonClasses(2).axonIds));
axonClasses(2).title = 'Thalamocortical axons';
axonClasses(2).tag = 'TC';
%}

% spine vs. shaft synapses
axonClasses(1).axonIds = find(conn.axonMeta.synCount);
axonClasses(1).synIds = find(synT.isSpine);
axonClasses(1).title = 'Spine synapses';
axonClasses(1).tag = 'Sp';

axonClasses(2).axonIds = find(conn.axonMeta.synCount);
axonClasses(2).synIds = find(~synT.isSpine & ~synT.isSoma);
axonClasses(2).title = 'Shaft synapses';
axonClasses(2).tag = 'Sh';

%% plot distribution of synapse size
binEdges = linspace(0, synAreaLim, 51);

fig = figure();
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
ax = axes(fig);
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
        'Normalization', 'probability', ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2);
end

ylabel(ax, 'Probability');
yticklabels(ax, arrayfun(@num2str, ...
    yticks(ax), 'UniformOutput', false));
xlabel(ax, 'Synapses per connection');

title( ...
    ax, {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
legend(ax, {axonClasses.title}, 'Location', 'NorthEast');

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

%% do comparisons again null hypothesis
controlCouplings = 2:5;
cvOf = @(d) std(d, 0, 2) ./ mean(d, 2);

for curCoupling = controlCouplings
    curFig = figure();
    curFig.Position(3:4) = [520, 756];
    
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

%% calculate baseline slope
synT = sortrows(synT, 'area', 'ascend');

rng(0);
randIds = floor(size(synT, 1) / 2);
randIds = randperm(size(synT, 1), 2 * randIds);

randIds = reshape(randIds, [], 2);
randSynAreas = synT.area(randIds);
randSynAreas = sort(randSynAreas, 2);

xLog = log10(randSynAreas(:, 2));
yLog = log10(randSynAreas(:, 1));

bl = [ones(numel(xLog), 1), xLog] \ yLog;
bl(1) = 10 ^ bl(1);

blFitF = @(x) bl(1) .* (x .^ bl(2));
blFitName = sprintf('Random pairs (y = %.2f x^{%.2f})', bl(1), bl(2));

%% look at doubly coupled neurites
[dupNeurites, ~, uniRows] = unique( ...
    synT(:, {'preAggloId', 'postAggloId'}), 'rows');

dupNeurites.synIds = accumarray( ...
    uniRows, synT.id, [], @(synIds) {transpose(synIds)});

uniCount = accumarray(uniRows, 1);
dupNeurites(uniCount ~= 2, :) = [];

% we now know that there are exactly two synapses
dupNeurites.synIds = cell2mat(dupNeurites.synIds);

% load synapse areas
[~, synRows] = ismember(dupNeurites.synIds, synT.id);
dupNeurites.synAreas = synT.area(synRows);

flipMask = ...
    dupNeurites.synAreas(:, 2) ...
  < dupNeurites.synAreas(:, 1);

dupNeurites.synIds(flipMask, :) = ...
    fliplr(dupNeurites.synIds(flipMask, :));
dupNeurites.synAreas(flipMask, :) = ...
    fliplr(dupNeurites.synAreas(flipMask, :));

[~, dupNeurites.interSynDist] = ismember( ...
    dupNeurites.preAggloId, interSyn.axonIds);
for curIdx = 1:size(dupNeurites, 1)
    curAxonIdx = dupNeurites.interSynDist(curIdx);
    curInterSynDists = interSyn.synToSynDists{curAxonIdx};
    
   [~, curSynIdx] = ismember( ...
        dupNeurites.synIds(curIdx, :), ...
        interSyn.synIds{curAxonIdx});
    dupNeurites.interSynDist(curIdx) = ...
        curInterSynDists(curSynIdx(1), curSynIdx(2));
end

%%
fig = figure();
ax = axes(fig);

hold(ax, 'on');
scatter(ax, ...
    dupNeurites.synAreas(:, 2), ...
    dupNeurites.synAreas(:, 1), 12, '+');

xlim([1E-2, 1E1]); xlabel('Axon-spine interface 1 (µm²)');
ylim([1E-2, 1E1]); ylabel('Axon-spine interface 2 (µm²)');

ax.XScale = 'log';
ax.YScale = 'log';

% do the fit
% taken from `matlab/+Analysis/+Script/bartolEtAl2015eLife.m` in `amotta`
xLog = log10(dupNeurites.synAreas(:, 2));
yLog = log10(dupNeurites.synAreas(:, 1));

b = [ones(numel(xLog), 1), xLog] \ yLog;
b(1) = 10 ^ b(1);

fitF = @(x) b(1) .* (x .^ b(2));
fitName = sprintf('Fit (y = %.2f x^{%.2f})', b(1), b(2));
rawName = sprintf('Raw data (n = %d)', numel(xLog));

fitRange = xlim();
fitRange = linspace(fitRange(1), fitRange(end), 2);
plot(fitRange, fitF(fitRange));
plot(fitRange, blFitF(fitRange));
plot(fitRange, fitRange, 'k--');

title( ...
   {'Same-axon same-dendrite spine synapses'; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
legend(rawName, fitName, blFitName, 'Location', 'NorthWest');

%% same-axon different-dendrite pairs
rng(0);

% get rid of any synapse size preference
saddT = synT(randperm(size(synT, 1)), :);

% TODO(amotta): If axon A makes multiple synapses onto dendrite D, we
% forget about all but one synapses for the following analysis. This makes
% the code easier. But I'm not yet sure whether this introduces some kind
% of bias...
[~, uniRows] = unique(saddT(:, {'preAggloId', 'postAggloId'}), 'rows');
saddT = saddT(uniRows, :);

% find axons that occur at least twice
[dendDupIds, ~, dendDupCount] = unique(saddT.preAggloId);
dendDupCount = accumarray(dendDupCount, 1);

dendDupIds(dendDupCount < 2) = [];
dendDupCount(dendDupCount < 2) = [];

assert(issorted(saddT.preAggloId));
[~, uniRows] = ismember(dendDupIds, saddT.preAggloId);

uniRows = arrayfun( ...
    @(r, c) r + (1:(2 * floor(c / 2)))' - 1, ...
	uniRows, dendDupCount, 'UniformOutput', false);
uniRows = cell2mat(uniRows);

saddT = saddT(uniRows, :);
saddT.pairId = ceil((1:size(saddT, 1)) / 2)';

% now that we've chosen the synapses, we can sort by size
saddT = sortrows(saddT, {'preAggloId', 'pairId', 'area'});

% sanity checks
assert(all(saddT.preAggloId(1:2:end) == saddT.preAggloId(2:2:end)));
assert(all(saddT.postAggloId(1:2:end) ~= saddT.postAggloId(2:2:end)));
assert(all(saddT.area(1:2:end) <= saddT.area(2:2:end)));

fig = figure();
ax = axes(fig);

hold(ax, 'on');
scatter(ax, ...
    saddT.area(2:2:end), ...
    saddT.area(1:2:end), 12, '+');

xlim([1E-2, 1E1]); xlabel('Axon-spine interface 1 (µm²)');
ylim([1E-2, 1E1]); ylabel('Axon-spine interface 2 (µm²)');

ax.XScale = 'log';
ax.YScale = 'log';

xLog = log10(saddT.area(2:2:end));
yLog = log10(saddT.area(1:2:end));

b = [ones(numel(xLog), 1), xLog] \ yLog;
b(1) = 10 ^ b(1);

fitF = @(x) b(1) .* (x .^ b(2));
fitName = sprintf('Fit (y = %.2f x^{%.2f})', b(1), b(2));
rawName = sprintf('Raw data (n = %d)', numel(xLog));

fitRange = xlim(ax);
fitRange = linspace(fitRange(1), fitRange(end), 2);
plot(ax, fitRange, fitF(fitRange));
plot(ax, fitRange, blFitF(fitRange));
plot(ax, fitRange, fitRange, 'k--');

title( ...
   {'Same-axon different-dendrite'; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
legend(rawName, fitName, blFitName, 'Location', 'NorthWest');

%% different-axon same-dendrite pairs
rng(0);

% get rid of any synapse size preference
dasdT = synT(randperm(size(synT, 1)), :);

% TODO(amotta): If axon A makes multiple synapses onto dendrite D, we
% forget about all but one synapses for the following analysis. This makes
% the code easier. But I'm not yet sure whether this introduces some kind
% of bias...
[~, uniRows] = unique(dasdT(:, {'postAggloId', 'preAggloId'}), 'rows');
dasdT = dasdT(uniRows, :);

% find axons that occur at least twice
[dendDupIds, ~, dendDupCount] = unique(dasdT.postAggloId);
dendDupCount = accumarray(dendDupCount, 1);

dendDupIds(dendDupCount < 2) = [];
dendDupCount(dendDupCount < 2) = [];

assert(issorted(dasdT.postAggloId));
[~, uniRows] = ismember(dendDupIds, dasdT.postAggloId);

uniRows = arrayfun( ...
    @(r, c) r + (1:(2 * floor(c / 2)))' - 1, ...
	uniRows, dendDupCount, 'UniformOutput', false);
uniRows = cell2mat(uniRows);

dasdT = dasdT(uniRows, :);
dasdT.pairId = ceil((1:size(dasdT, 1)) / 2)';

% now that we've chosen the synapses, we can sort by size
dasdT = sortrows(dasdT, {'postAggloId', 'pairId', 'area'});

% sanity checks
assert(all(dasdT.postAggloId(1:2:end) == dasdT.postAggloId(2:2:end)));
assert(all(dasdT.preAggloId(1:2:end) ~= dasdT.preAggloId(2:2:end)));
assert(all(dasdT.area(1:2:end) <= dasdT.area(2:2:end)));

fig = figure();
ax = axes(fig);

hold(ax, 'on');
scatter(ax, ...
    dasdT.area(2:2:end), ...
    dasdT.area(1:2:end), 12, '+');

xlim([1E-2, 1E1]); xlabel('Axon-spine interface 1 (µm²)');
ylim([1E-2, 1E1]); ylabel('Axon-spine interface 2 (µm²)');

ax.XScale = 'log';
ax.YScale = 'log';

xLog = log10(dasdT.area(2:2:end));
yLog = log10(dasdT.area(1:2:end));

b = [ones(numel(xLog), 1), xLog] \ yLog;
b(1) = 10 ^ b(1);

fitF = @(x) b(1) .* (x .^ b(2));
fitName = sprintf('Fit (y = %.2f x^{%.2f})', b(1), b(2));
rawName = sprintf('Raw data (n = %d)', numel(xLog));

fitRange = xlim(ax);
fitRange = linspace(fitRange(1), fitRange(end), 2);
plot(ax, fitRange, fitF(fitRange));
plot(ax, fitRange, blFitF(fitRange));
plot(ax, fitRange, fitRange, 'k--');

title( ...
   {'Different-axon same-dendrite'; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
legend(rawName, fitName, blFitName, 'Location', 'NorthWest');

%% different-axon different-dendrite
rng(0);

% get rid of any synapse size preference
daddT = synT(randperm(size(synT, 1)), :);

[~, uniRows] = unique(daddT(:, {'preAggloId', 'postAggloId'}), 'rows');
daddT = daddT(uniRows, :);

curPos = 0;
curSynsOpen = true(size(daddT, 1), 1);
curIds = zeros(0, 2);

while curPos < height(daddT)
    curPos = curPos + 1;
    
    if ~curSynsOpen(curPos)
        continue;
    end
    
    curSyn = daddT(curPos, :);
    curPossIds = all(bsxfun(@ne, ...
        [daddT.preAggloId, daddT.postAggloId], ...
        [curSyn.preAggloId, curSyn.postAggloId]), 2);
    curPossIds = find(curPossIds & curSynsOpen);
    
    if ~isempty(curPossIds)
        curPossId = curPossIds(randi(numel(curPossIds)));
        curSynsOpen(cat(2, curPossId, curPos)) = false;
        curIds(end + 1, :) = [curPos, curPossId]; %#ok
    end
end

curIds = reshape(transpose(curIds), [], 1);
daddT = daddT(curIds, :);

daddT.pairId = ceil((1:size(daddT, 1)) / 2)';
daddT = sortrows(daddT, {'pairId', 'area'});

% sanity checks
assert(all(daddT.postAggloId(1:2:end) ~= daddT.postAggloId(2:2:end)));
assert(all(daddT.preAggloId(1:2:end) ~= daddT.preAggloId(2:2:end)));
assert(all(daddT.area(1:2:end) <= daddT.area(2:2:end)));

fig = figure();
ax = axes(fig);

hold(ax, 'on');
scatter(ax, ...
    daddT.area(2:2:end), ...
    daddT.area(1:2:end), 12, '+');

xlim([1E-2, 1E1]); xlabel('Axon-spine interface 1 (µm²)');
ylim([1E-2, 1E1]); ylabel('Axon-spine interface 2 (µm²)');

ax.XScale = 'log';
ax.YScale = 'log';

xLog = log10(daddT.area(2:2:end));
yLog = log10(daddT.area(1:2:end));

b = [ones(numel(xLog), 1), xLog] \ yLog;
b(1) = 10 ^ b(1);

fitF = @(x) b(1) .* (x .^ b(2));
fitName = sprintf('Fit (y = %.2f x^{%.2f})', b(1), b(2));
rawName = sprintf('Raw data (n = %d)', numel(xLog));

fitRange = xlim(ax);
fitRange = linspace(fitRange(1), fitRange(end), 2);
plot(ax, fitRange, fitF(fitRange));
plot(ax, fitRange, blFitF(fitRange));
plot(ax, fitRange, fitRange, 'k--');

title( ...
   {'Different-axon different-dendrite'; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
legend(rawName, fitName, blFitName, 'Location', 'NorthWest');

%%
maxCv = 1.5;

cvData = struct;
cvData(1).name = 'Random pairs';
cvData(1).data = randSynAreas;

cvData(2).name = 'Same-axon same-dendrite';
cvData(2).data = dupNeurites.synAreas;

cvData(3).name = 'Same-axon different-dendrite';
cvData(3).data = reshape(saddT.area, 2, [])';

cvData(4).name = 'Different-axon same-dendrite';
cvData(4).data = reshape(dasdT.area, 2, [])';

cvData(5).name = 'Different-axon different-dendrite';
cvData(5).data = reshape(daddT.area, 2, [])';

fig = figure();
ax = axes(fig);
hold(ax, 'on');

for curIdx = 1:numel(cvData)
    cvData(curIdx).cv = ...
        std(cvData(curIdx).data, 0, 2) ...
        ./ mean(cvData(curIdx).data, 2);
    
    cvData(curIdx).median = ...
        median(cvData(curIdx).cv);
    
   [~, cvData(curIdx).pValue] = ...
        kstest2(cvData(curIdx).cv, cvData(1).cv);
    
    if curIdx == 1; continue; end
    
    cvData(curIdx).name = sprintf( ...
        '%s (p = %.2g)', ...
        cvData(curIdx).name, ...
        cvData(curIdx).pValue);
end

for curIdx = 1:numel(cvData)
    histogram( ...
        ax, cvData(curIdx).cv, linspace(0, maxCv, 51), ...
        'Normalization', 'probability', 'DisplayStyle', 'stairs');
end

title( ...
   {'Synapse size variability'; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
legend(ax, cvData.name, 'Location', 'SouthWest');

ax.TickDir = 'out';
ax.XLim = [0, maxCv];
xlabel(ax, 'Coefficient of variation');
ylabel(ax, 'Probability');

%% estimate circuit learnedness
% TODO(amotta): justify + automatically derive
threshCv = 0.66;

areaB = mean(cvData(2).cv > threshCv);
areaC = mean(cvData(1).cv < threshCv);
areaA = mean(cvData(2).cv < threshCv) - areaC;

% NOTE(amotta): areas A + B + C = 1
learnLowerBound = areaA / (areaA + areaB + areaC) %#ok
learnUpperBound = (1 - areaB) / (areaA + areaB + areaC) %#ok

%% plot same-axon same-dendrite CV vs. intersynapse distance
maxInterSynDistUm = 25;

% calculate inter-synapse distance
dupNeurites.synAreaCv = ...
    std(dupNeurites.synAreas, 0, 2) ...
    ./ mean(dupNeurites.synAreas, 2);

fig = figure();
ax = axes(fig);
hold(ax, 'on');

scatter(ax, ...
    dupNeurites.interSynDist / 1E3, ...
    dupNeurites.synAreaCv, 'x');

mask = (dupNeurites.interSynDist < maxInterSynDistUm * 1E3);
b = [ ...
    ones(size(dupNeurites.interSynDist(mask))), ...
    dupNeurites.interSynDist(mask) / 1E3] \ dupNeurites.synAreaCv(mask);

fitF = @(x) b(1) + x .* b(2);
rawName = sprintf('Raw data (n = %d)', numel(xLog));
fitName = sprintf( ...
    'Linear fit (y = %.2g + %.4gd; d_{max} = %g µm)', ...
    b(1), b(2), maxInterSynDistUm);

fitRange = xlim(ax);
fitRange = linspace(fitRange(1), fitRange(end), 2);
plot(ax, fitRange, fitF(fitRange));

ax.TickDir = 'out';
xlabel(ax, 'Intersynapse distance (µm)');
ylabel(ax, 'Coefficient of variation');
legend(ax, rawName, fitName, 'Location', 'SouthEast');

title( ...
   {'Distance dependence of synaptic consistency'; ...
    info.git_repos{1}.hash}, 'FontWeight', 'normal', 'FontSize', 10);
