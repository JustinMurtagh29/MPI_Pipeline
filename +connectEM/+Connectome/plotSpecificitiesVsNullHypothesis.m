% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

targetClasses = { ...
    'Somata', 'ProximalDendrite', 'ApicalDendrite', ...
    'SmoothDendrite', 'AxonInitialSegment', 'OtherDendrite'};

minSynPre = 10;

info = Util.runInfo();
Util.showRunInfo(info);

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, ~, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);

[conn, axonClasses] = ...
    connectEM.Connectome.prepareForSpecificityAnalysis( ...
        conn, axonClasses, 'minSynPre', minSynPre);

classConnectome = ...
    connectEM.Connectome.buildClassConnectome( ...
        conn, 'targetClasses', targetClasses);

%% calculate target class innervation probabilities for null model
clear cur*;

for curIdx = 1:numel(axonClasses)
    curNullProbs = classConnectome(axonClasses(curIdx).nullAxonIds, :);
    curNullProbs = connectEM.Specificity.calcFirstHitProbs(curNullProbs, 'multinomial');
    axonClasses(curIdx).nullTargetClassProbs = curNullProbs;
end

%% show first hit probabilities
firstHitProbs = cat(1, axonClasses.nullTargetClassProbs);
firstHitProbs = array2table(firstHitProbs, 'VariableNames', targetClasses);
firstHitProbs.Properties.RowNames = {axonClasses.title};

fprintf('# First-hit probabilities\n');
disp(firstHitProbs);

%% plot
clear cur*;
curAxonClasses = axonClasses;

% NOTE(amotta): The list of axons to induce the null model is no longer
% needed. Let's remove it to cause an error in case a weird code path still
% tries to use it.
curAxonClasses = curAxonClasses(1:4);
curAxonClasses = rmfield(curAxonClasses, 'nullAxonIds');

for curIdx = 1:numel(curAxonClasses)
    plotAxonClass( ...
        info, classConnectome, ...
        targetClasses, curAxonClasses(curIdx));
end

%% plotting
function plotAxonClass(info, classConn, targetClasses, axonClass)
    axonSpecs = classConn(axonClass.axonIds, :);
    synCounts = sum(axonSpecs, 2);
    axonSpecs = axonSpecs ./ synCounts;
    
    %% preparations
    nullTargetClassProbs = axonClass.nullTargetClassProbs;
    axonNullProbs = connectEM.Specificity.calcChanceProbs( ...
        classConn, axonClass.axonIds, nullTargetClassProbs, ...
        'distribution', 'binomial');
    
    %% plotting
    fig = figure;
    fig.Color = 'white';
    fig.Position(3:4) = [1850, 1025];
    
    binEdges = linspace(0, 1, 21);
    axes = cell(size(targetClasses));
    pValAxes = cell(size(targetClasses));

    for classIdx = 1:numel(targetClasses)
        className = targetClasses{classIdx};
        classProb = nullTargetClassProbs(classIdx);
        
        axonClassSpecs = axonSpecs(:, classIdx);
        axonClassNullProbs = axonNullProbs(:, classIdx);
        
        % Null hypothesis
       [nullSynFrac, nullAxonCount] = ...
            connectEM.Specificity.calcExpectedDist( ...
                synCounts, classProb, 'distribution', 'binomial');
            
        ksProb = ...
            connectEM.Specificity.kolmogorovSmirnovTest( ...
                axonClassSpecs, nullSynFrac, ...
                'nullWeights', nullAxonCount, ...
                'tail', 'smaller');
        
        nullBinId = discretize(nullSynFrac, binEdges);
        nullBinCount = accumarray(nullBinId, nullAxonCount);

        % Measured
        ax = subplot(3, numel(targetClasses), classIdx);
        axis(ax, 'square');
        hold(ax, 'on');
        
        histogram(ax, ...
            axonClassSpecs, ...
            'BinEdges', binEdges, ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2, ...
            'FaceAlpha', 1);
        histogram(ax, ...
            'BinEdges', binEdges, ...
            'BinCounts', nullBinCount, ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2, ...
            'FaceAlpha', 1);

        xlabel(ax, 'Synapse fraction');
        ax.XAxis.TickDirection = 'out';
        ax.XAxis.Limits = [0, 1];
        
        ylabel(ax, 'Axons');
        ax.YAxis.TickDirection = 'out';
        ax.YAxis.Limits(1) = 10 ^ (-0.1);
        ax.YAxis.Scale = 'log';
        
        title(ax, ...
            {className; sprintf('p = %g (tailed KS)', ksProb)}, ...
            'FontWeight', 'normal', 'FontSize', 10);
        axes{classIdx} = ax;
        
        %% p-values
        curBinEdges = linspace(-1E-3, 1 + 1E-3, numel(binEdges));
        
       [expChanceProbs, expChanceCounts] = ...
            connectEM.Specificity.calcExpectedChanceProbDist( ...
                synCounts, classProb);
            
        curExpCounts = discretize(expChanceProbs, curBinEdges);
        curExpCounts = accumarray(curExpCounts, expChanceCounts);
            
        ax = subplot( ...
            3, numel(targetClasses), ...
            numel(targetClasses) + classIdx);
        axis(ax, 'square');
        hold(ax, 'on');
        
        histogram(ax, ...
            axonClassNullProbs, ...
            'BinEdges', curBinEdges, ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2, ...
            'FaceAlpha', 1);
        histogram(ax, ...
            'BinCounts', curExpCounts, ...
            'BinEdges', curBinEdges, ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2, ...
            'FaceAlpha', 1);
        
        ax.XLim = curBinEdges([1, end]);
        
        pValAxes{classIdx} = ax;
        
        %% alternative visualization
        % Compare p-value distribution against expectation:
        % We'd expect there to be `theta` percent of axons with a p-value
        % below `theta`. If there are, however, significantly more axons
        % with a p-value below `theta`, something interesting is going on.
        curPVal = sort(axonClassNullProbs, 'ascend');
        curPVal = reshape(curPVal, 1, []);
        
        ax = subplot( ...
            3, numel(targetClasses), ...
            2 * numel(targetClasses) + classIdx);
        axis(ax, 'square');
        hold(ax, 'on');
        
       [curPVal, ~, curPAxonFrac] = unique(curPVal);
        curPAxonFrac = accumarray(curPAxonFrac, 1);
        curPAxonFrac = cumsum(curPAxonFrac) / sum(curPAxonFrac);
        
        % Conservative estimate of false detection rate (FDR)
        curFdrEst = cumsum(expChanceCounts);
        curFdrEst = curFdrEst / curFdrEst(end);
        
        curFdrEst = interp1(expChanceProbs, curFdrEst, curPVal);
        curFdrEst = curFdrEst(:) ./ curPAxonFrac(:);
        
        curThetaIdx = 1 + find( ...
            curFdrEst(1:(end - 1)) <= 0.2 ...
          & curFdrEst(2:end) > 0.2, 1);
        
        plot(ax, curPVal, curFdrEst, 'LineWidth', 1);
        
        xlim(ax, [0, 1]);
        ylim(ax, [0, 1.2]);
        xlabel(ax, 'p-value');
        ylabel(ax, 'Estimated FDR');
        
        if ~isempty(curThetaIdx)
            curThetaPVal = curPVal(curThetaIdx);
            
            plot(ax, ...
                repelem(curThetaPVal, 2), ylim(ax), ...
                'Color', 'black', 'LineStyle', '--');
            title(ax, ...
                sprintf('p = %f', curThetaPVal), ...
                'FontWeight', 'normal', 'FontSize', 10);
        end
    end
    
    % Legend
    ax = axes{end};
    axPos = ax.Position;
    leg = legend(ax, ...
        'Observed', ...
        'Binomial model', ...
        'Location', 'East');
    leg.Box = 'off';
    
    % Fix positions
    ax.Position = axPos;
    leg.Position(1) = sum(axPos([1, 3])) + 0.005;

    axes = horzcat(axes{:});
    yMax = max(arrayfun(@(a) a.YAxis.Limits(end), axes));
    for ax = axes; ax.YAxis.Limits(end) = yMax; end
    
    annotation( ...
        'textbox', [0, 0.9, 1, 0.1], ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
        'String', { ...
            'Observed synapse fractions vs. null hypothesis'; ...
            axonClass.title; info.git_repos{1}.hash});
end
