% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
%   Benedikt Staffler <benedikt.staffler@brain.mpg.de>
%
% Differences to _v1:
%   - synapses: SynapseAgglos_v6_somaH.mat and the corresponding connectome

clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', ...
    'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v6-somaH.mat');

targetClasses = { ...
    'Somata', 'ProximalDendrite', 'ApicalDendrite', ...
    'SmoothDendrite', 'AxonInitialSegment', 'OtherDendrite'};

minSynPre = 10;
info = Util.runInfo();

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, ~, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);

% Inhibitory whole cell → smooth dendrite
% Excitatory whole cell → proximal dendrite
wcMask = conn.denMeta.targetClass == 'WholeCell';
inMask = conn.denMeta.isInterneuron;

conn.denMeta.targetClass(wcMask &  inMask) = 'SmoothDendrite';
conn.denMeta.targetClass(wcMask & ~inMask) = 'ProximalDendrite';

classConnectome = ...
    connectEM.Connectome.buildClassConnectome( ...
        conn, 'targetClasses', targetClasses);
    
%% generate a class with all axons
allAxonClass = struct;
allAxonClass.axonIds = find( ...
    conn.axonMeta.synCount >= minSynPre);
allAxonClass.nullAxonIds = find( ...
    conn.axonMeta.synCount >= minSynPre);
allAxonClass.title = sprintf( ...
    'all axons with ≥ %d synapses (n = %d)', ...
    minSynPre, numel(allAxonClass.axonIds));

axonClasses(end + 1) = allAxonClass;

%% plot
for curIdx = 1:numel(axonClasses)
    outlier = plotAxonClass( ...
        info, conn.axonMeta, classConnectome, ...
        targetClasses, axonClasses(curIdx));
end


%% plotting
function outlier = plotAxonClass(info, axonMeta, classConn, targetClasses, axonClass)
    axonSpecs = classConn(axonClass.axonIds, :);
    axonSpecs = axonSpecs ./ sum(axonSpecs, 2);
    
    
    %% preparations
    
    %     nullModel = 'binomial';
    nullModel = 'drmn';
    
    % calculate overall synapse probabilities
    targetClassSyns = sum(classConn(axonClass.nullAxonIds, :), 1);
    targetClassProbs = targetClassSyns / sum(targetClassSyns);
    
    % dirichlet-multinomial fit (polya)
    try
        a = polya_fit_simple(classConn(axonClass.axonIds, :)); % dirichlet-multinomial dist params % requires fastfit toolbox https://github.com/tminka/fastfit
        fitted_polya = true;
    catch err
        fitted_polya = false;
    end
    
    if fitted_polya
        % polya fit samples (maximum likelihood)
        pol_s = polya_sample(a, sum(classConn(axonClass.axonIds, :), 2));
        pol_s = pol_s ./ sum(pol_s, 2);
    %     pol_m =  sum(targetClassSyns) .* a ./ sum(a);

        % multinomial samples (maximum likelihood)
        mn_s = mnrnd(sum(classConn(axonClass.axonIds, :), 2), targetClassProbs);
        mn_s = mn_s ./ sum(mn_s, 2);
        
%         % fit for axons that actually do innervate a target class
%         figure;
%         i = 1;
%         histogram(axonSpecs(axonSpecs(:,i) > 0, 1), -0.025:0.05:1.25)
%         hold on
%         histogram(pol_s(pol_s(:,i) > 0, 1), -0.025:0.05:1.25)
%         histogram(mn_s(mn_s(:,i) > 0, 1), -0.025:0.05:1.25)
%         legend('Observed', 'Polya', 'Multinomial');
%         title(sprintf('Synapse fraction histogram (target %d - innervations only)', i))
    end
    
    switch nullModel
        case 'binomial'
            axonNullProbs = connectEM.Specificity.calcChanceProbs( ...
                classConn, axonClass.axonIds, axonClass.nullAxonIds, ...
                'distribution', 'binomial');
        case 'drmn'
            axonNullProbs = connectEM.Specificity.calcChanceProbs( ...
                classConn, axonClass.axonIds, axonClass.nullAxonIds, ...
                'distribution', 'drmn', 'alpha', a);
    end
    
    %% plotting
    fig = figure;
    fig.Color = 'white';
    fig.Position(3:4) = [1850, 885];
    
    binEdges = linspace(0, 1, 21);
    axes = cell(size(targetClasses));
    pValAxes = cell(size(targetClasses));
    
    outlier = cell(length(targetClasses), 1);

    for classIdx = 1:numel(targetClasses)
        className = targetClasses{classIdx};
        classProb = targetClassProbs(classIdx);
        
        axonClassSpecs = axonSpecs(:, classIdx);
        axonClassNullProbs = axonNullProbs(:, classIdx);
        
        % Null hypothesis
       [nullSynFrac, nullAxonCount] = ...
            connectEM.Specificity.calcExpectedDist( ...
                axonMeta.synCount(axonClass.axonIds), ...
                classProb, 'distribution', 'binomial');
            
        nullBinId = discretize(nullSynFrac, binEdges);
        nullBinCount = accumarray(nullBinId, nullAxonCount);
            
        % null hypothesis dirichlet multinomial marginal (beta-binomial)
        if fitted_polya
            p1 = a(classIdx);
            p2 = sum(a(setdiff(1:length(a), classIdx)));
            
            warning('off', 'all')
            [nullSynFrac, nullAxonCount] = ...
                connectEM.Specificity.calcExpectedDist( ...
                    axonMeta.synCount(axonClass.axonIds), ...
                    classProb, 'distribution', 'bb', 'a', p1, 'b', p2);
            warning('on', 'all')
            
            nullBinId_bb = discretize(nullSynFrac, binEdges);
            nullBinCount_bb = accumarray(nullBinId_bb, nullAxonCount);
        end

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
        
        if fitted_polya
            histogram(ax, ...
                'BinEdges', binEdges, ...
                'BinCounts', nullBinCount_bb, ...
                'DisplayStyle', 'stairs', ...
                'LineWidth', 2, ...
                'FaceAlpha', 1);
        end

        xlabel(ax, 'Synapse fraction');
        ax.XAxis.TickDirection = 'out';
        ax.XAxis.Limits = [0, 1];
        
        ylabel(ax, 'Axons');
        ax.YAxis.TickDirection = 'out';
        ax.YAxis.Limits(1) = 10 ^ (-0.1);
        ax.YAxis.Scale = 'log';
        
        title(ax, className, 'FontWeight', 'normal', 'FontSize', 10);
        axes{classIdx} = ax;
        
        %% p-values
        ax = subplot( ...
            3, numel(targetClasses), ...
            numel(targetClasses) + classIdx);
        axis(ax, 'square');
        hold(ax, 'on');
        
        % p-value distribution under null hypothesis (binomial)
        [p_vals_h0, ax_count_h0] = ...
            connectEM.Specificity.calcExpectedChanceProbDist( ...
            axonMeta.synCount(axonClass.axonIds), classProb);
        ax_count_h0 = cumsum(ax_count_h0);
        ax_count_h0 = ax_count_h0 ./ ax_count_h0(end);
        
        if fitted_polya
            % p-value distribution under null hypothesis (dirichlet-multinomial)
            p1 = a(classIdx);
            p2 = sum(a(setdiff(1:length(a), classIdx)));
            pdf = @(n)Math.Prob.bbinopdf(0:n, n, p1, p2);
            warning('off', 'all');
            [p_vals_h0_bb, ax_count_h0_bb] = ...
                connectEM.Specificity.calcExpectedChanceProbDist( ...
                axonMeta.synCount(axonClass.axonIds), classProb, 'pdf', pdf);
            warning('on', 'all');
            ax_count_h0_bb = cumsum(ax_count_h0_bb);
            ax_count_h0_bb = ax_count_h0_bb ./ ax_count_h0_bb(end);
            
            % get outliers defined as axons for which the upper cdf is
            % smaller than some threshold out_t according to the null model
            % which is equivalent to saying that the p-value is smaller
            % than out_t
            out_t = 0.001;
            outlier{classIdx} = find(axonClassNullProbs < out_t);
        end
        
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
        hold(ax, 'on');
        plot( ...
            ax, binEdges([1, end]), [1, 1], ...
            'Color', ax.ColorOrder(2, :));
        pl = plot( ...
            ax, curPVal, curPRatio, ...
            'Color', ax.ColorOrder(1, :), ...
            'LineWidth', 2);
        
        if ~isempty(curThetaIdx)
            plot( ...
                ax, curPVal([curThetaIdx, curThetaIdx]), ax.YLim, ...
                'Color', 'black', 'LineStyle', '--');
            title(ax, ...
                sprintf('p = %.2f', curPVal(curThetaIdx)), ...
                'FontWeight', 'normal', 'FontSize', 10);
        end
        
        ax.TickDir = 'out';
        xlabel(ax, 'p-value');
        xlim(ax, binEdges([1, end]));
        ylim(ax, [0, 2]);
        
        pValAxes{classIdx} = ax;
        
        if classIdx == numel(targetClasses)
            legend(pl, ...
                sprintf('observed p-value under %s', nullModel), ...
                'Location', 'NorthEast');
        end
        
        %% alternative visualization
        ax = subplot( ...
            3, numel(targetClasses), ...
            2 * numel(targetClasses) + classIdx);
        axis(ax, 'square');
        hold(ax, 'on');
        
        curPAxonFrac = linspace(0, 1, numel(curPVal));
        
        plot(curPVal, curPAxonFrac);
        pl = plot(p_vals_h0, ax_count_h0);
        if fitted_polya
            pl2 = plot(p_vals_h0_bb, ax_count_h0_bb);
            pl = [pl, pl2];
        end
%         plot([0, 1], [0, 1]);
        
        if ~isempty(curThetaIdx)
            plot(ax, ...
                curPVal([curThetaIdx, curThetaIdx]), ax.YLim, ...
                'Color', 'black', 'LineStyle', '--');
        end
        
        ax.XLim = [0, 1];
        xlabel(ax, 'p-value');
        ylabel(ax, {'Fraction of axons'; 'with p < x'});
        
        if classIdx == numel(targetClasses)
            legend(pl, ...
                'sampled p-value under binomial', ...
                'sampled p-value under drmn', ...
                'Location', 'NorthEast');
        end
    end
    
    % Legend
    ax = axes{end};
    axPos = ax.Position;
    leg = legend(ax, ...
        'Observed', ...
        'Binomial model', ...
        'Location', 'East');
    if fitted_polya
        leg = legend(ax, ...
            'Observed', ...
            'Binomial model', ...
            'Dirichlet-Multinomial', ...
            'Location', 'East');
    end
    leg.Box = 'off';
    
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
