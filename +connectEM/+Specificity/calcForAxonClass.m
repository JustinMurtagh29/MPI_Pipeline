function specs = calcForAxonClass( ...
        classConn, targetClasses, axonClass, varargin)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opts = struct;
    opts.info = [];
    opts.showPlot = false;
    opts.fdrThresh = 0.2;
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
            curFdrEst(1:(end - 1)) <= opts.fdrThresh ...
          & curFdrEst(2:end) > opts.fdrThresh, 1);
      
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
