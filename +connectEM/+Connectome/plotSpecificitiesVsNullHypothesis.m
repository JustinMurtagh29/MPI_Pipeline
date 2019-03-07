% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

targetClasses = { ...
    'Somata', 'ProximalDendrite', 'ApicalDendrite', ...
    'SmoothDendrite', 'AxonInitialSegment', 'OtherDendrite'};

% NOTE(amotta): Rows correspond to shaft and spine synapses, respectively.
% Columns correspond to the target classes listed above. All-zeros mark
% rates that were not measured (because biologically implausible).
inhFpRates = [ ...
     (1 / 20), (1 / 19), (5 / 17), (1 / 17), (4 / 18), (3 / 17);  % shaft
     (000000), (6 / 18), (6 / 18), (4 / 16), (000000), (6 / 17)]; % spine

minSynPre = 10;

info = Util.runInfo();
Util.showRunInfo(info);

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, syn, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);
[axonClasses.tag] = deal('Exc', 'Inh', 'TC', 'CC', 'Unclear');

%% build class connectome for shaft synapses
% This is needed to simulate the effect of FP inhibitory synapses.
clear cur*;

curShaftConn = conn;
curShaftConn.connectome.synIdx = cellfun( ...
    @(ids) ids(syn.synapses.type(ids) == 'Shaft'), ...
    curShaftConn.connectome.synIdx, 'UniformOutput', false);

curShaftConn = ...
    connectEM.Connectome.prepareForSpecificityAnalysis( ...
        curShaftConn, axonClasses, 'minSynPre', minSynPre);
classConnShafts = ...
    connectEM.Connectome.buildClassConnectome( ...
        curShaftConn, 'targetClasses', targetClasses);
    
%% build complete class connectome
clear cur*;

[conn, axonClasses] = ...
    connectEM.Connectome.prepareForSpecificityAnalysis( ...
        conn, axonClasses, 'minSynPre', minSynPre);

classConn = ...
    connectEM.Connectome.buildClassConnectome( ...
        conn, 'targetClasses', targetClasses);
    
%% prepare for analysis while accounting for false positive inh. synapses
%{
[axonClasses.synFalsePosRates] = deal([]);
[axonClasses.synFalsePosMethod] = deal([]);

axonClasses(2).synFalsePosRates = inhFpRates;
axonClasses(2).synFalsePosMethod = 'binomial';
%}

%% plot
clear cur*;

evalT = cell(size(axonClasses));
barFracs = cell(size(axonClasses));
pValThreshs = [0.025, 0.05, 0.1, 0.2, 0.3];

for curClassId = 1:numel(axonClasses)
    curAxonClass = axonClasses(curClassId);
    curFpRates = [];
    
    %{
    curFpRates = curAxonClass.synFalsePosRates;
    curFpMethod = curAxonClass.synFalsePosMethod;
    
    % NOTE(amotta): If empirically determined synapse false positive rates
    % have been specified, we simulate the effect of synapse removal under
    % a binomial distribution. This necessitates multiple multiple runs to
    % get a feeling for the variability across trials.
    switch isempty(curFpRates)
        case true, numRuns = 1;
        case false, numRuns = 20;
    end
    %}
    
    numRuns = numel(pValThreshs);
    
    curEvalT = array2table( ...
        nan(numRuns, 1 + numel(targetClasses)), ...
        'VariableNames', [{'Overall'}, targetClasses]);
    curBarFracs = nan(numRuns, 2 * numel(targetClasses) - 1);
    
    for curRun = 1:numRuns
        curPValThresh = pValThreshs(curRun);
        curConn = classConn(curAxonClass.nullAxonIds, :);
        
        if ~isempty(curFpRates)
            % Take into account false positive synapse rates
            curShaftConn = classConnShafts(curAxonClass.nullAxonIds, :);
            curSpineConn = curConn - curShaftConn;

            rng(curRun);
            
            switch curFpMethod
              case 'binomial'
                % Binomial correction for false positive synapses
                curShaftConn = curShaftConn - binornd(curShaftConn, ...
                    repelem(curFpRates(1, :), size(curShaftConn, 1), 1));
                curSpineConn = curSpineConn - binornd(curSpineConn, ...
                    repelem(curFpRates(2, :), size(curSpineConn, 1), 1));
            
              case 'constFraction'
                % NOTE(amotta): Correction for false positive synapses
                % by removal of a random subset of synapses. With the
                % binomial approach the number of remove synapses
                % changes across runs. This method here keeps the
                % fraction of removed synapses constant.
                curShaftConn = ...
                    removeConstantFractionOfSynapses( ...
                        curShaftConn, curFpRates(1, :));
                curSpineConn = ...
                    removeConstantFractionOfSynapses( ...
                        curSpineConn, curFpRates(2, :));

              otherwise
                error('Invalid method "%s"', curFpMethod);
            end
            
            curConn = curShaftConn + curSpineConn;
        end
        
        curAxonClass.nullTargetClassProbs = ...
            connectEM.Specificity.calcFirstHitProbs( ...
                curConn, 'oneVersusRestBinomial');
        curAxonClass.specs = plotAxonClass( ...
            classConn, targetClasses, curAxonClass, ...
            'pValThresh', curPValThresh, 'showPlot', false, 'info', info);
        
       [curA, curB, curC] = ...
            connectEM.Specificity.calcAxonFractions( ...
                curAxonClass, targetClasses);
        
        curEvalT{curRun, :} = [curA, curB];
        curBarFracs(curRun, :) = curC;
    end
    
    evalT{curClassId} = curEvalT;
    barFracs{curClassId} = curBarFracs;
end

%% quantitative evaluation
clear cur*;

for curIdx = 1:numel(axonClasses)
    curAxonClass = axonClasses(curIdx);
    
    curEvalT = evalT{curIdx};
    curEvalT.pValThresh = pValThreshs(:);
    curEvalT = circshift(curEvalT, 1, 2);
    
    fprintf('%s\n\n', curAxonClass.title);
    disp(curEvalT); fprintf('\n');
end

%% Plot fraction of axons specific
% Copy & paste from +connectEM/+Figure/fractionOfAxonsSpecific.m
clear cur*;
plotClasses = 1:2;

curFig = figure();
curAx = axes(curFig);
hold(curAx, 'on');

for curIdx = 1:numel(plotClasses)
    curId = plotClasses(curIdx);
    curFracs = evalT{curId}.Overall;
    
    plot(curAx, repelem(curIdx, numel(curFracs)), curFracs');
end

curLines = curAx.Children;
set(curLines, 'Color', 'black', 'Marker', '.', 'MarkerSize', 16);

curAx.XLim = [0.5, numel(plotClasses) + 0.5];
curAx.YLim = [0, 1];

curAx.XTick = 1:numel(plotClasses);
curAx.XTickLabel = {axonClasses(plotClasses).tag};
ylabel(curAx, {'Fraction of'; 'axons specific'});

connectEM.Figure.config(curFig, info);
curFig.Position(3:4) = [125, 160];

%% plotting
function specs = plotAxonClass( ...
        classConn, targetClasses, axonClass, varargin)
    opts = struct;
    opts.info = [];
    opts.showPlot = false;
    opts.pValThresh = 0.2;
    opts = Util.modifyStruct(opts, varargin{:});
    
    specs = struct;
    
    nullTargetClassProbs = axonClass.nullTargetClassProbs;
    axonSpecs = classConn(axonClass.axonIds, :);
    synCounts = sum(axonSpecs, 2);
    axonSpecs = axonSpecs ./ synCounts;
    
    %% preparations
    axonNullProbs = ...
        connectEM.Specificity.calcChanceProbs( ...
            classConn, axonClass.axonIds, nullTargetClassProbs, ...
            'distribution', 'binomial');
    
    binEdges = linspace(0, 1, 21);
    
    if opts.showPlot
        fig = figure;
        axes = cell(size(targetClasses));
        pValAxes = cell(size(targetClasses));
    end
    
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
        
        % p-values
        curBinEdges = linspace(-1E-3, 1 + 1E-3, numel(binEdges));
        
       [expChanceProbs, expChanceCounts] = ...
            connectEM.Specificity.calcExpectedChanceProbDist( ...
                synCounts, classProb);
        
        curExpCounts = discretize(expChanceProbs, curBinEdges);
        curExpCounts = accumarray(curExpCounts, expChanceCounts);
        
        % alternative visualization
        % Compare p-value distribution against expectation:
        % We'd expect there to be `theta` percent of axons with a p-value
        % below `theta`. If there are, however, significantly more axons
        % with a p-value below `theta`, something interesting is going on.
        curPVal = sort(axonClassNullProbs, 'ascend');
        curPVal = reshape(curPVal, 1, []);
        
       [curPVal, ~, curPAxonFrac] = unique(curPVal);
        curPAxonFrac = accumarray(curPAxonFrac, 1);
        curPAxonFrac = cumsum(curPAxonFrac) / sum(curPAxonFrac);
        
        % Conservative estimate of false detection rate (FDR)
        curFdrEst = cumsum(expChanceCounts);
        curFdrEst = curFdrEst / curFdrEst(end);
        
        curFdrEst = interp1(expChanceProbs, curFdrEst, curPVal);
        curFdrEst = curFdrEst(:) ./ curPAxonFrac(:);
        
        curThetaPVal = 1 + find( ...
            curFdrEst(1:(end - 1)) <= opts.pValThresh ...
          & curFdrEst(2:end) > opts.pValThresh, 1);
      
        curThetaPVal = curPVal(curThetaPVal);
        if isempty(curThetaPVal); curThetaPVal = nan; end
        
        if ksProb < 0.01 ...
                && not(isnan(curThetaPVal)) ...
                && not(strcmpi(className, 'OtherDendrite'))
            curSpecs = struct;
            curSpecs.pThresh = curThetaPVal;
            curSpecs.axonIds = axonClass.axonIds( ...
                axonClassNullProbs < curThetaPVal);
            specs.(className) = curSpecs;
        end

        %% Plotting
        if ~opts.showPlot; continue; end
        
        ax = subplot(3, numel(targetClasses), classIdx);
        axis(ax, 'square');
        hold(ax, 'on');
        
        histogram(ax, ...
            axonClassSpecs, ...
            'BinEdges', binEdges, ...
            'DisplayStyle', 'stairs');
        histogram(ax, ...
            'BinEdges', binEdges, ...
            'BinCounts', nullBinCount, ...
            'DisplayStyle', 'stairs');

        xlabel(ax, 'Synapse fraction');
        ax.XAxis.TickDirection = 'out';
        ax.XAxis.Limits = [0, 1];
        
        ylabel(ax, 'Axons');
        ax.YAxis.TickDirection = 'out';
        ax.YAxis.Limits(1) = 10 ^ (-0.1);
        ax.YAxis.Scale = 'log';
        
        title(ax, ...
            {className; sprintf('p = %g (tailed KS)', ksProb)});
        axes{classIdx} = ax;
            
        ax = subplot( ...
            3, numel(targetClasses), ...
            numel(targetClasses) + classIdx);
        axis(ax, 'square');
        hold(ax, 'on');
        
        histogram(ax, ...
            axonClassNullProbs, ...
            'BinEdges', curBinEdges, ...
            'DisplayStyle', 'stairs');
        histogram(ax, ...
            'BinCounts', curExpCounts, ...
            'BinEdges', curBinEdges, ...
            'DisplayStyle', 'stairs');
        
        pValAxes{classIdx} = ax;
        
        ax = subplot( ...
            3, numel(targetClasses), ...
            2 * numel(targetClasses) + classIdx);
        axis(ax, 'square');
        hold(ax, 'on');
        
        plot(ax, curPVal, curFdrEst, 'LineWidth', 1);
        
        xlim(ax, [0, 1]);
        ylim(ax, [0, 1.2]);
        xlabel(ax, 'p-value');
        ylabel(ax, 'Estimated FDR');
        
        if ~isnan(curThetaPVal)
            plot(ax, ...
                repelem(curThetaPVal, 2), ylim(ax), ...
                'Color', 'black', 'LineStyle', '--');
            title(ax, sprintf('p = %f', curThetaPVal));
        end
    end
    
    if ~opts.showPlot; return; end
    
    % Legend
    ax = axes{end};
    axPos = ax.Position;
    leg = legend(ax, ...
        'Observed', ...
        'Binomial model', ...
        'Location', 'East');
    
    % Fix positions
    ax.Position = axPos;
    leg.Position(1) = sum(axPos([1, 3])) + 0.005;

    axes = horzcat(axes{:});
    yMax = max(arrayfun(@(a) a.YAxis.Limits(end), axes));
    for ax = axes; ax.YAxis.Limits(end) = yMax; end
    
    curTitle = { ...
        'Observed synapse fractions vs. null hypothesis'; axonClass.title};
    annotation(fig, 'textbox', [0, 0.9, 1, 0.1], 'String', curTitle);
    
    connectEM.Figure.config(fig, opts.info);
    fig.Position(3:4) = [1850, 1150];
end

function classConn = ...
        removeConstantFractionOfSynapses(classConn, removeRate)
    assert(all(0 <= removeRate(:) & removeRate(:) <= 1));
    classCount = size(classConn, 2);
    
    if isscalar(removeRate)
        removeRate = repelem(removeRate, 1, classCount);
    else
        % NOTE(amotta): This will (intentionally) throw an error if the
        % number of elements in `removeRate` does not agree with the number
        % of target classes.
        removeRate = reshape(removeRate, 1, classCount);
    end
    
    for curIdx = 1:classCount
        curRemoveRate = removeRate(curIdx);
        
        % Find non-empty entries in class connectome
       [curIds, ~, curCounts] = find(classConn(:, curIdx));
        curIds = repelem(curIds(:), curCounts(:));
        
        % Select subset of synapses to drop
        curRemoveCount = round(curRemoveRate * numel(curIds));
        curIds = curIds(randperm(numel(curIds), curRemoveCount));
        
        % Update class connectome
       [curIds, ~, curCounts] = unique(curIds);
        curCounts = accumarray(curCounts, 1);
        
        classConn(curIds, curIdx) = ...
            classConn(curIds, curIdx) - curCounts;
    end
end
