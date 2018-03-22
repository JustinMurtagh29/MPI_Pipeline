% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
param = struct;
param.saveFolder = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connName = 'connectome_ax18a_deWC01wSp';

% For recovery of all TC axon definitions
connOldName = 'connectome_axons_18_a_ax_spine_syn_clust';
eqClassesFile = fullfile(param.saveFolder, 'SVGDB', 'agglos', 'ax18a_deWC01wSp', 'eClass.mat');

minSynPre = 10;
info = Util.runInfo();

%% loading data
% Load all TC axon definitions
connOld = connectEM.Connectome.load(param, connOldName);
tcAxonsOld = connOld.axons(connOld.axonMeta.axonClass == 'Thalamocortical');

eqClasses = load(eqClassesFile);
eqClasses = eqClasses.idsToEClass;

% TODO(amotta): Talk to Benedikt and find out why `eqClasses(segIds)` does
% not consist of a single non-zero number.
tcAxonsNew = cellfun( ...
    @(ids) mode(nonzeros(eqClasses(ids))), ...
    tcAxonsOld, 'UniformOutput', false);

assert(all(cellfun(@isscalar, tcAxonsNew)));
tcAxonsNew = cell2mat(tcAxonsNew);

conn = connectEM.Connectome.load(param, connName);
[classConnectome, targetClasses] = ...
    connectEM.Connectome.buildClassConnectome(conn);
axonClasses = ...
    connectEM.Connectome.buildAxonClasses(conn, 'minSynPre', minSynPre);

% Overwrite definition of TC axons
axonClasses(3).axonIds = tcAxonsNew;
axonClasses(3).title = sprintf( ...
   ['thalamocortical axons (n = %d; ', ...
    'from old connectome)'], numel(tcAxonsNew));

%% plot
for curIdx = 1:numel(axonClasses)
    plotAxonClass( ...
        info, conn.axonMeta, classConnectome, ...
        targetClasses, axonClasses(curIdx));
end

%% plotting
function plotAxonClass(info, axonMeta, classConn, targetClasses, axonClass)
    axonSpecs = classConn(axonClass.axonIds, :);
    axonSpecs = axonSpecs ./ sum(axonSpecs, 2);
    
    %% preparations
    axonNullProbs = connectEM.Specificity.calcChanceProbs( ...
        classConn, axonClass.axonIds, axonClass.nullAxonIds, ...
        'distribution', 'binomial');
    
    % calculate overall synapse probabilities
    targetClassSyns = sum(classConn(axonClass.nullAxonIds, :), 1);
    targetClassProbs = targetClassSyns / sum(targetClassSyns);
    
    %% plotting
    fig = figure;
    fig.Color = 'white';
    fig.Position(3:4) = [1850, 800];
    
    binEdges = linspace(0, 1, 21);
    axes = cell(size(targetClasses));
    pValAxes = cell(size(targetClasses));

    for classIdx = 1:numel(targetClasses)
        className = char(targetClasses(classIdx));
        classProb = targetClassProbs(classIdx);
        
        axonClassSpecs = axonSpecs(:, classIdx);
        axonClassNullProbs = axonNullProbs(:, classIdx);
        isSpecific = axonClassNullProbs < 0.01;
        
        % Null hypothesis
       [nullSynFrac, nullAxonCount] = ...
            connectEM.Specificity.calcExpectedDist( ...
                axonMeta.synCount(axonClass.axonIds), ...
                classProb, 'distribution', 'binomial');
        
        nullBinId = discretize(nullSynFrac, binEdges);
        nullBinCount = accumarray(nullBinId, nullAxonCount);

        % Measured
        ax = subplot(2, numel(targetClasses), classIdx);
        axis(ax, 'square');
        hold(ax, 'on');

        histogram(ax, ...
            axonClassSpecs(isSpecific), ...
            'BinEdges', binEdges, ...
            'EdgeColor', 'none', ...
            'FaceAlpha', 1);
        histogram(ax, ...
            axonClassSpecs, ...
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
        
        % Show number of overly specific axons. Two approaches are used:
        % 1. Find a synapse fraction threshold above which less than one
        %    axon is expected under the null hypothesis. Then count the
        %    number of axons above this threshold.
        % 2. Count the number of axons which are above the null hypothesis
        %    distribution. This does not truly represent overly-specific
        %    axons.
        overBinId = 1 + find(nullBinCount > 1, 1, 'last');
        
        obsBinCount = histcounts(axonClassSpecs, binEdges);
        overThreshCount = sum(obsBinCount(overBinId:end));
        specificCount = sum(isSpecific);
        
        title(ax, { ...
            sprintf('%d with S ≥ %.1f', ...
                overThreshCount, binEdges(overBinId)); ...
            sprintf('%d with p ≤ 1 %%', specificCount)}, ...
            'FontWeight', 'normal', 'FontSize', 10);
        
        axes{classIdx} = ax;
        
        %% p-values
        ax = subplot( ...
            2, numel(targetClasses), ...
            numel(targetClasses) + classIdx);
        axis(ax, 'square');
        hold(ax, 'on');
        
        % Compare p-value distribution against expectation:
        % We'd expect there to be `theta` percent of axons with a p-value
        % below `theta`. If there are, however, significantly more axons
        % with a p-value below `theta`, something interesting is going on.
        curPVal = sort(axonClassNullProbs, 'ascend');
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
            axonClassNullProbs, binEdges, ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2);
        xlim(ax, binEdges([1, end]));
        ylim(ax, [0, size(axonSpecs, 1)]);
        
        xlabel(ax, 'p-value');
        ax.XAxis.TickDirection = 'out';
       [ax.YAxis.TickDirection] = deal('out');
        
        pValAxes{classIdx} = ax;
    end
    
    % Legend
    ax = axes{end};
    axPos = ax.Position;
    leg = legend(ax, ...
        'p < 1 %', ...
        'Observed', ...
        'Binomial model', ...
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
            axonClass.title; info.git_repos{1}.hash});
end