% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a_ax_spine_syn_clust.mat');
outDir = '/home/amotta/Desktop';

showPlots = true;
info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

points = Seg.Global.getSegToPointMap(param);

conn = load(connFile);

%% Based on spine synapse fraction
% This is justified according to Yunfeng's data
minSynCount = 10;
spineSynFracThresh = 0.5;

candIds = find( ...
    conn.denMeta.synCount >= minSynCount ...
  & conn.denMeta.targetClass ~= 'Somata');

spineSynFrac = ...
    conn.denMeta.spineSynCount ...
    ./ conn.denMeta.synCount;

if showPlots
    binEdges = linspace(0, 1, 51);
    
    fig = figure();
    fig.Color = 'white';
    fig.Position(3:4) = [600, 600];
    
    ax = axes(fig);
    axis(ax, 'square');
    
    yyaxis(ax, 'right');
    histogram(ax, ...
        spineSynFrac(candIds), binEdges, ...
        'Normalization', 'cdf', ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2);
    ylabel(ax, 'Cumulative distribution function');
    
    yyaxis(ax, 'left');
    hold(ax, 'on');
    histogram(ax, ...
        spineSynFrac(candIds), binEdges, ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2);
    plot(ax, ...
        repelem(spineSynFracThresh, 1, 2), ...
        ylim(ax), 'Color', 'black', 'LineStyle', '--');
    ylabel(ax, 'Dendrites');
    
    ax.TickDir = 'out';
    xlim(ax, binEdges([1, end]));
    xlabel(ax, 'Spine synapse fraction');
    
    annotation(fig, ...
        'textbox', [0, 0.9, 1, 0.1], ...
        'String', { ...
            info.filename; info.git_repos{1}.hash;
            sprintf( ...
                '%d dendrites with >= %d synapses', ...
                numel(candIds), minSynCount)}, ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center');
end

smoothIds = candIds;
smoothIds(spineSynFrac(smoothIds) > spineSynFracThresh) = [];

%% Export to webKNOSSOS
rng(0);
randIds = randperm(numel(smoothIds));
randIds = smoothIds(randIds(1:25));
randIds = reshape(randIds, 1, []);

dendPoints = cellfun( ...
    @(segIds) points(segIds, :), ...
    conn.dendrites(randIds), ...
    'UniformOutput', false);

dendNames = arrayfun( ...
    @(idx, id) sprintf( ...
        '%0*d. Agglomerate %d', ...
        ceil(log10(1 + numel(randIds))), idx, id), ...
	1:numel(randIds), randIds, 'UniformOutput', false);

skel = Skeleton.fromMST(dendPoints, param.raw.voxelSize);
skel.names = dendNames;

skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));
skel.write(fullfile(outDir, 'smooth-dendrite-candidates.nml'));