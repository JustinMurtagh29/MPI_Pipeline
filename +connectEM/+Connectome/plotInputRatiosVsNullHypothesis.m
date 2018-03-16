% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
param = struct;
param.saveFolder = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connName = 'connectome_axons_18_a_ax_spine_syn_clust';

minSynPost = 10;
info = Util.runInfo();

%% loading data
conn = ...
    connectEM.Connectome.load(param, connName);
axonClasses = unique(conn.axonMeta.axonClass);

%% build class connectome
classConnectome = ...
    connectEM.Connectome.buildClassConnectome( ...
        conn, 'targetClasses', [], 'axonClasses', axonClasses);

classConnectome = transpose(classConnectome);
axonClasses = reshape(axonClasses, 1, []);

%% build dendrite class(es)
dendClasses = struct;
dendClasses(1).ids = find( ...
    conn.denMeta.targetClass ~= 'Somata' ...
  & conn.denMeta.targetClass ~= 'AxonInitialSegment' ...
  & conn.denMeta.synCount >= minSynPost);
dendClasses(1).nullIds = dendClasses(1).ids;
dendClasses(1).title = sprintf( ...
    'Dendrites with ≥ %d synapses (n = %d)', ...
    minSynPost, numel(dendClasses(1).ids));
    
%% plot
for curIdx = 1:numel(dendClasses)
    plotAxonClass( ...
        info, classConnectome, ...
        axonClasses, dendClasses(curIdx));
end

%% plotting
function plotAxonClass(info, classConn, axonClasses, dendClass)
    import connectEM.Specificity.calcExpectedFractionDist;
    classConn = classConn(dendClass.nullIds, :);
    
    % Plotting
    fig = figure();
    fig.Color = 'white';
    fig.Position(3:4) = [840, 440];
    binEdges = linspace(0, 1, 21);
    
    % Plot expected e / (e + i)
    excId = find(axonClasses == 'Excitatory');
    inhId = find(axonClasses == 'Inhibitory');
    
    obsFrac = classConn(:, [excId, inhId]);
    obsFrac = obsFrac(:, 1) ./ sum(obsFrac, 2);
    obsFrac(isnan(obsFrac)) = 0;
    
   [expFrac, expCount] = ...
        calcExpectedFractionDist(classConn, excId, inhId);
        
	binId = discretize(expFrac, binEdges);
    binCount = accumarray(binId, expCount);
    
    ax = subplot(1, 2, 1);
    axis(ax, 'square');
    hold(ax, 'on');
    
    histogram(ax, ...
        obsFrac, ...
        'BinEdges', binEdges, ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2);
    histogram(ax, ...
        'BinEdges', binEdges, ...
        'BinCounts', binCount, ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2);
    
    xlabel(ax, 'exc / (exc + inh)');
    xlim(ax, binEdges([1, end]));
    
    ax.YAxis.Limits(1) = 10 ^ (-0.1);
    ax.YScale = 'log';
    ax.TickDir = 'out';
    
    % Plot expected tc / (tc + e)
    tcId = find(axonClasses == 'Thalamocortical');
    excId = find(axonClasses == 'Excitatory');
    
    obsFrac = classConn(:, [tcId, excId]);
    obsFrac = obsFrac(:, 1) ./ sum(obsFrac, 2);
    obsFrac(isnan(obsFrac)) = 0;
    
   [expFrac, expCount] = ...
        calcExpectedFractionDist(classConn, tcId, excId);
        
	binId = discretize(expFrac, binEdges);
    binCount = accumarray(binId, expCount);
    
    ax = subplot(1, 2, 2);
    axis(ax, 'square');
    hold(ax, 'on');
    
    histogram(ax, ...
        obsFrac, ...
        'BinEdges', binEdges, ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2);
    histogram(ax, ...
        'BinEdges', binEdges, ...
        'BinCounts', binCount, ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2);
    
    xlabel(ax, 'tc / (tc + exc)');
    xlim(ax, binEdges([1, end]));
    
    ax.YAxis.Limits(1) = 10 ^ (-0.1);
    ax.YScale = 'log';
    ax.TickDir = 'out';
    
    % Legend
    axPos = ax.Position;
    leg = legend(ax, ...
        'Observed', ...
        'Expected (multinomial)', ...
        'Location', 'NorthEast');
    ax.Position = axPos;

    annotation(fig, ...
        'textbox', [0, 0.9, 1, 0.1], ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'center', ...
        'String', { ...
            info.filename; ...
            info.git_repos{1}.hash; ...
            dendClass.title});
end