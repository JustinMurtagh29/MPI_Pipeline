% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');

% Control with linearized axons
ctrlConnFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_b_linearized_ax_spine_syn_clust.mat');
ctrlSynFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v3_ax_spine_clustered.mat');

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, syn, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);

[~, synToSynFile] = fileparts(connFile);
synToSynFile = sprintf('%s_synToSynDists.mat', synToSynFile);
synToSynFile = fullfile(fileparts(connFile), synToSynFile);
synToSyn = load(synToSynFile);

ctrlConn = load(ctrlConnFile);
ctrlSyn = load(ctrlSynFile);

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

connectEM.Consistency.plotSizeHistogram( ...
    info, synT, curPlotConfigs, ...
    'binEdges', linspace(-2, +1, 61), ...
    'scale', 'log10');

%% Illustrate synapse size similarity
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
for curPlotConfig = plotConfigs
    curPairConfigs = ...
        connectEM.Consistency.buildPairConfigs(synT, curPlotConfig);
    curFig = ...
        connectEM.Consistency.plotVariabilityHistogram( ...
            info, synT, curPlotConfig, curPairConfigs(:));
    curFig.Position(3:4) = [370, 540];
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
curMaxDistUm = 20;

curPlotConfig = plotConfigs(1);
curPairConfig = connectEM.Consistency.buildPairConfigs(synT, curPlotConfig);
curPairConfig = curPairConfig(1);

curT = struct2table(rmfield(curPairConfig, 'title'));
curT.preAggloId = synT.preAggloId(curT.synIdPairs(:, 1));
curT.postAggloId = synT.postAggloId(curT.synIdPairs(:, 1));

curT.preDist = zeros(size(curT.preAggloId));
curT.postDist = zeros(size(curT.postAggloId));
for curIdx = 1:numel(curT.preDist)
    curPreAggloId = curT.preAggloId(curIdx);
    curPostAggloId = curT.postAggloId(curIdx);
    curSynIds = synT.id(curT.synIdPairs(curIdx, :));
    
    % Distance along axonal side
    curPreSynIds = synToSyn.axonSynIds{curPreAggloId};
   [~, curPreSynIds] = ismember(curSynIds, curPreSynIds);
    curPreDist = synToSyn.axonSynToSynDists{curPreAggloId};
    curPreDist = curPreDist(curPreSynIds(1), curPreSynIds(2));
    
    % Distance along axonal side
    curPostSynIds = synToSyn.dendSynIds{curPostAggloId};
   [~, curPostSynIds] = ismember(curSynIds, curPostSynIds);
    curPostDist = synToSyn.dendSynToSynDists{curPostAggloId};
    curPostDist = curPostDist(curPostSynIds(1), curPostSynIds(2));
    
    curT.preDist(curIdx) = curPreDist;
    curT.postDist(curIdx) = curPostDist;
end

curT.cv = synT.area(curT.synIdPairs);
curT.cv = std(curT.cv, 0, 2) ./ mean(curT.cv, 2);

curPreT = curT(curT.preDist / 1E3 < curMaxDistUm, :);
curFitPre = fit(curPreT.preDist / 1E3, curPreT.cv, 'poly1');

curPostT = curT(curT.postDist / 1E3 < curMaxDistUm, :);
curFitPost = fit(curPostT.postDist / 1E3, curPostT.cv, 'poly1');

% Plot
fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [560, 360];

% Presynaptic side
ax = subplot(1, 2, 1);
axis(ax, 'square');
hold(ax, 'on');

plot(ax, curPreT.preDist / 1E3, curPreT.cv, '.');
plot(ax, ax.XLim, curFitPre(ax.XLim), 'Color', 'black');
title(ax, {'Presynaptic', sprintf( ...
    'CV = %.2f + %.4fd', curFitPre.p2, curFitPre.p1)}, ...
    'FontWeight', 'normal', 'FontSize', 10);

% Postsynaptic side
ax = subplot(1, 2, 2);
axis(ax, 'square');
hold(ax, 'on');

plot(ax, curPostT.postDist / 1E3, curPostT.cv, '.');
plot(ax, ax.XLim, curFitPost(ax.XLim), 'Color', 'black');
title(ax, {'Postsynaptic', sprintf( ...
    'CV = %.2f + %.4fd', curFitPost.p2, curFitPost.p1)}, ...
    'FontWeight', 'normal', 'FontSize', 10);

axes = fig.Children;
[ax.TickDir] = deal('out');

ax = axes(end);
xlabel(ax, 'Intersynapse distance (µm)');
ylabel(ax, 'Coefficient of variation');

annotation(fig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'String', { ...
        info.filename; ...
        info.git_repos{1}.hash; ...
        curPairConfig.title}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');

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
