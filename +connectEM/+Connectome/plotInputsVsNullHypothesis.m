% NOTE
%   This script is based on `plotSpecificitiesVsNullHypothesis.m`.
%
%   If you find a bug in this script, please check whether the above file
%   suffers from the same problem! (Yes, I know that this screams "pack it
%   into a resusable function already.")
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');

wcFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto-and-manual.mat');

splitNmlDir = fileparts(fileparts(mfilename('fullpath')));
splitNmlDir = fullfile(splitNmlDir, '+WholeCell', '+Script', 'annotations');

minSynPost = 10;
info = Util.runInfo();

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, syn] = connectEM.Connectome.load(param, connFile);
conn = connectEM.Connectome.prepareForSpecificityAnalysis(conn);
axonClasses = unique(conn.axonMeta.axonClass);

wcData = load(wcFile);
splitNmlT = connectEM.WholeCell.loadSplitNmls(splitNmlDir);

%% split whole cells
dendT = table;
dendT.cellId = repelem( ...
    splitNmlT.aggloId, cellfun(@numel, splitNmlT.dendNodes));
dendT.dendId = cell2mat(cellfun( ...
    @(dendNodes) transpose(1:numel(dendNodes)), ...
    splitNmlT.dendNodes, 'UniformOutput', false));

dendT.agglo = wcData.dendrites(dendT.cellId);
dendT.cellId = wcData.idxWholeCells(dendT.cellId);
dendT.agglo = arrayfun( ...
    @(a, nC) unique(rmmissing(a.nodes(nC{1}, 4))), ...
    dendT.agglo, cat(1, splitNmlT.dendNodes{:}), ...
    'UniformOutput', false);

synIds = repelem( ...
    transpose(1:height(syn.synapses)), ...
    cellfun(@numel, syn.synapses.postsynId));
synSegIds = cell2mat(syn.synapses.postsynId);

dendT.classConn = cellfun( ...
    @(segIds) unique(synIds(ismember(synSegIds, segIds))), ...
    dendT.agglo, 'UniformOutput', false);
clear synIds synSegIds;

synAxonIds = repelem( ...
    conn.connectome.edges(:, 1), ...
    cellfun(@numel, conn.connectome.synIdx));
synAxonIds = double(conn.axonMeta.axonClass(synAxonIds));
synIds = cell2mat(conn.connectome.synIdx);

synAxonIds = cell2mat(accumarray( ...
    synIds, synAxonIds, [size(syn.synapses, 1), 1], ...
    @(axonIds) {accumarray(axonIds, 1, [4, 1])'}, {zeros(1, 4)}));
dendT.classConn = cell2mat(cellfun( ...
    @(ids) sum(synAxonIds(ids, :), 1), ...
    dendT.classConn, 'UniformOutput', false));
clear synIds synAxonIds;

% Get rid of dendrites with less than 50 synapses
dendT(sum(dendT.classConn, 2) < 50, :) = [];

%% build class connectome
classConnectome = ...
    connectEM.Connectome.buildClassConnectome( ...
        conn, 'targetClasses', [], 'axonClasses', axonClasses);

[dendMeta, classConnectome] = ...
    connectEM.Connectome.prepareForFullCellInputAnalysis( ...
        conn.denMeta, classConnectome);

classConnectome = transpose(classConnectome);
axonClasses = reshape(axonClasses, 1, []);

%% build dendrite class(es)
dendClasses = struct;
dendClasses(1).ids = find( ...
    dendMeta.targetClass ~= 'Somata' ...
  & dendMeta.targetClass ~= 'AxonInitialSegment' ...
  & dendMeta.targetClass ~= 'FullInput' ...
  & dendMeta.synCount >= minSynPost);
dendClasses(1).nullIds = dendClasses(1).ids;
dendClasses(1).title = sprintf( ...
    'Dendrites with ≥ %d synapses (n = %d)', ...
    minSynPost, numel(dendClasses(1).ids));

dendClasses(2).ids = find( ...
    dendMeta.targetClass == 'FullInput' ...
  & dendMeta.synCount >= 500);
dendClasses(2).nullIds = dendClasses(2).ids;
dendClasses(2).title = sprintf( ...
    'Whole cells with ≥ %d synapses (n = %d)', ...
    500, numel(dendClasses(2).ids));

%% plot
for curIdx = 1:numel(dendClasses)
    plotAxonClass(info, classConnectome, axonClasses, dendClasses(curIdx));
end

%%
curDendClass = struct;
curDendClass.ids = transpose(1:height(dendT));
curDendClass.nullIds = curDendClass.ids;
curDendClass.title = sprintf( ...
    'Dendrites with ≥ 50 synapses (n = %d)', ...
    numel(curDendClass.ids));
    
plotAxonClass(info, dendT.classConn, axonClasses, curDendClass);

%% plot variability across dendrites of a cell
%{
for curId = reshape(unique(dendT.cellId), 1, [])
    curDendT = dendT(dendT.cellId == curId, :);
    if height(curDendT) < 2; continue; end
    
    curDendClass = struct;
    curDendClass.ids = transpose(1:height(curDendT));
    curDendClass.nullIds = curDendClass.ids;
    curDendClass.title = sprintf( ...
        'Dendrites of cell %d (n = %d)', ...
        curId, height(curDendT));
    
    plotAxonClass(info, curDendT.classConn, axonClasses, curDendClass);
end
%}

%% plotting
function plotAxonClass(info, classConn, axonClasses, dendClass)
    dendSynFracs = classConn(dendClass.ids, :);
    dendSynFracs = dendSynFracs ./ sum(dendSynFracs, 2);
    
    %% preparations
    nullProbs = connectEM.Specificity.calcChanceProbs( ...
        classConn, dendClass.ids, dendClass.nullIds, ...
        'distribution', 'binomial');
    
    % calculate overall synapse probabilities
    axonClassProbs = sum(classConn(dendClass.nullIds, :), 1);
    axonClassProbs = axonClassProbs / sum(axonClassProbs);
    
    %% plotting
    fig = figure;
    fig.Color = 'white';
    fig.Position(3:4) = [1850, 900];
    
    binEdges = linspace(0, 1, 51);
    axes = cell(size(axonClasses));
    pValAxes = cell(size(axonClasses));

    for classIdx = 1:numel(axonClasses)
        className = char(axonClasses(classIdx));
        axonClassProb = axonClassProbs(classIdx);
        
        dendClassSynFracs = dendSynFracs(:, classIdx);
        dendClassNullProbs = nullProbs(:, classIdx);
        
        % Null hypothesis
       [nullSynFrac, nullDendCount] = ...
            connectEM.Specificity.calcExpectedDist( ...
                sum(classConn(dendClass.ids, :), 2), ...
                axonClassProb, 'distribution', 'binomial');
        
        nullBinId = discretize(nullSynFrac, binEdges);
        nullBinCount = accumarray(nullBinId, nullDendCount);
        
        % Measured
        ax = subplot(2, numel(axonClasses), classIdx);
        axis(ax, 'square');
        hold(ax, 'on');
        histogram(ax, ...
            dendClassSynFracs, ...
            'BinEdges', binEdges, ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2);
        histogram(ax, ...
            'BinEdges', binEdges, ...
            'BinCounts', nullBinCount, ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2);
        
        xlabel(ax, className);
        ax.XAxis.TickDirection = 'out';
        ax.XAxis.Limits = [0, 1];
        
        %{
        ax.YAxis.TickDirection = 'out';
        ax.YAxis.Limits(1) = 10 ^ (-0.1);
        ax.YAxis.Scale = 'log';
        %}
        axes{classIdx} = ax;
        
        %% p-values
        ax = subplot( ...
            2, numel(axonClasses), ...
            numel(axonClasses) + classIdx);
        axis(ax, 'square');
        hold(ax, 'on');
        
        % Compare p-value distribution against expectation:
        % We'd expect there to be `theta` percent of dendrites with a
        % p-value below `theta`. If there are, however, significantly more
        % dendrites with a p-value below `theta`, something interesting is
        % going on.
        curPVal = sort(dendClassNullProbs, 'ascend');
        curPVal = reshape(curPVal, 1, []);

        curPRatio = (1:numel(curPVal)) ./ numel(curPVal);
        curPRatio = curPRatio ./ curPVal;
        
        % Find chance level (i.e., ratio 1)
        curThetaIdx = find(curPRatio > 1, 1);
        
        % No threshold if chance level was reached from below.
        if mean(curPRatio(1:curThetaIdx) > 1) < 0.5
            curThetaIdx = [];
        end
        
        if ~isempty(curThetaIdx)
            curThetaIdx = curThetaIdx - 1 + find( ...
                curPRatio(curThetaIdx:end) < 1, 1);
        end
        
        % Plotting
        yyaxis(ax, 'right');
        hold(ax, 'on');
        plot(ax, curPVal, curPRatio, 'LineWidth', 2);
        plot(ax, binEdges([1, end]), [1, 1], 'LineStyle', '--');
        ylim(ax, [0, 2]);
        
        if ~isempty(curThetaIdx)
            plot(ax, ...
                curPVal([curThetaIdx, curThetaIdx]), ...
                ax.YLim, 'LineStyle', '--');
            title(ax, ...
                sprintf('p = %.2f', curPVal(curThetaIdx)), ...
                'FontWeight', 'normal', 'FontSize', 10);
        end
        
        yyaxis(ax, 'left');
        histogram(ax, ...
            dendClassNullProbs, binEdges, ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2);
        xlim(ax, binEdges([1, end]));
        ylim(ax, [0, size(dendSynFracs, 1)]);
        
        xlabel(ax, 'p-value');
        ax.XAxis.TickDirection = 'out';
       [ax.YAxis.TickDirection] = deal('out');
        
        pValAxes{classIdx} = ax;
    end
    
    % Legend
    ax = axes{end};
    axPos = ax.Position;
    leg = legend(ax, ...
        'Observed', ...
        'Expected (binomial)', ...
        'Location', 'East');
    
    % Fix positions
    ax.Position = axPos;
    leg.Position(1) = sum(axPos([1, 3])) + 0.005;

    axes = horzcat(axes{:});
    yMax = max(arrayfun(@(a) a.YAxis.Limits(end), axes));
    for ax = axes; ax.YAxis.Limits(end) = yMax; end
    
    pValAxes = horzcat(pValAxes{:});
    yMax = max(arrayfun(@(a) a.YAxis(1).Limits(end), pValAxes));
   [pValAxes.YLim] = deal([0, yMax]);

    annotation( ...
        'textbox', [0, 0.9, 1, 0.1], ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
        'String', { ...
            'Observed synapse fractions vs. null hypothesis'; ...
            dendClass.title; info.git_repos{1}.hash});
end
