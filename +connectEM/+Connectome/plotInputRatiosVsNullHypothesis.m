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

%{
dendClasses(2).ids = find( ...
    conn.denMeta.targetClass == 'Somata' ...
  & conn.denMeta.synCount >= minSynPost);
dendClasses(2).nullIds = dendClasses(2).ids;
dendClasses(2).title = sprintf( ...
    'Somata with ≥ %d synapses (n = %d)', ...
    minSynPost, numel(dendClasses(2).ids));

dendClasses(3).ids = find( ...
    conn.denMeta.targetClass == 'WholeCell' ...
  & conn.denMeta.synCount >= minSynPost);
dendClasses(3).nullIds = dendClasses(3).ids;
dendClasses(3).title = sprintf( ...
    'Whole cells with ≥ %d synapses (n = %d)', ...
    minSynPost, numel(dendClasses(3).ids));
%}
    
%% plot
for curIdx = 1:numel(dendClasses)
    plotAxonClass( ...
        info, classConnectome, ...
        axonClasses, dendClasses(curIdx));
end

%% plotting
function plotAxonClass(info, classConn, axonClasses, dendClass)
    fig = figure();
    binEdges = linspace(0, 1, 21);
    
    % Plot expected e / (e + i)
    excId = find(axonClasses == 'Excitatory');
    inhId = find(axonClasses == 'Inhibitory');
    
   [expFrac, expCount] = ...
        connectEM.Specificity.calcExpectedFractionDist( ...
            classConn(dendClass.nullIds, :), excId, inhId); %#ok
        
	binId = discretize(expFrac, binEdges);
    binCount = accumarray(binId, expCount);
    
    ax = subplot(1, 2, 1);
    histogram(ax, ...
        'BinEdges', binEdges, ...
        'BinCounts', binCount, ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2);
    
    xlabel(ax, 'exc / (exc + inh)');
    xlim(ax, binEdges([1, end]));
    ax.TickDir = 'out';
    
    % Plot expected tc / (tc + e)
    tcId = find(axonClasses == 'Thalamocortical');
    excId = find(axonClasses == 'Excitatory');
    
   [expFrac, expCount] = ...
        connectEM.Specificity.calcExpectedFractionDist( ...
            classConn(dendClass.nullIds, :), tcId, excId); %#ok
        
	binId = discretize(expFrac, binEdges);
    binCount = accumarray(binId, expCount);
    
    ax = subplot(1, 2, 2);
    histogram(ax, ...
        'BinEdges', binEdges, ...
        'BinCounts', binCount, ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2);
    
    xlabel(ax, 'tc / (tc + exc)');
    xlim(ax, binEdges([1, end]));
    ax.TickDir = 'out';
end