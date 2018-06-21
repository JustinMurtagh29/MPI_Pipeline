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

minSynPost = 10;
info = Util.runInfo();

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

conn = connectEM.Connectome.load(param, connFile);
conn = connectEM.Connectome.prepareForSpecificityAnalysis(conn);
axonClasses = unique(conn.axonMeta.axonClass);

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
    plotAxonClass( ...
        info, dendMeta, classConnectome, ...
        axonClasses, dendClasses(curIdx));
end

%% plotting
function plotAxonClass(info, dendMeta, classConn, targetClasses, dendClass)
    dendSynFracs = classConn(dendClass.ids, :);
    dendSynFracs = dendSynFracs ./ sum(dendSynFracs, 2);
    
    % calculate overall synapse probabilities
    axonClassProbs = classConn(dendClass.nullIds, :);
    axonClassProbs = axonClassProbs / sum(axonClassProbs, 2);
    axonClassProbs = mean(axonClassProbs, 1);
    
    %% preparations
    nullProbs = connectEM.Specificity.calcChanceProbs( ...
        classConn, dendClass.ids, axonClassProbs, ...
        'distribution', 'binomial');
    
    %% plotting
    fig = figure;
    fig.Color = 'white';
    fig.Position(3:4) = [1850, 900];
    
    binEdges = linspace(0, 1, 51);
    axes = cell(size(targetClasses));
    pValAxes = cell(size(targetClasses));

    for classIdx = 1:numel(targetClasses)
        className = char(targetClasses(classIdx));
        axonClassProb = axonClassProbs(classIdx);
        
        dendClassSynFracs = dendSynFracs(:, classIdx);
        dendClassNullProbs = nullProbs(:, classIdx);
        
        % Null hypothesis
       [nullSynFrac, nullDendCount] = ...
            connectEM.Specificity.calcExpectedDist( ...
                dendMeta.synCount(dendClass.ids), ...
                axonClassProb, 'distribution', 'binomial');
        
        nullBinId = discretize(nullSynFrac, binEdges);
        nullBinCount = accumarray(nullBinId, nullDendCount);
        
        % Measured
        ax = subplot(2, numel(targetClasses), classIdx);
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

        ax.YAxis.TickDirection = 'out';
        ax.YAxis.Limits(1) = 10 ^ (-0.1);
        ax.YAxis.Scale = 'log';
        
        axes{classIdx} = ax;
        
        %% p-values
        ax = subplot( ...
            2, numel(targetClasses), ...
            numel(targetClasses) + classIdx);
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
