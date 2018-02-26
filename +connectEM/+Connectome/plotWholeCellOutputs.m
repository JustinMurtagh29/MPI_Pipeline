% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
outputDir = '/home/amotta/Desktop/whole-cell-output-distributions';

wcFile = fullfile(rootDir, 'aggloState', 'wholeCells_08.mat');
synFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v3_ax_spine_clustered.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a_ax_spine_syn_clust.mat');

info = Util.runInfo();

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);
segMass = Seg.Global.getSegToSizeMap(param);
segCentroids = Seg.Global.getSegToCentroidMap(param);

syn = load(synFile);

wc = load(wcFile);
wc = wc.wholeCells;

[~, connName] = fileparts(connFile);
conn = connectEM.Connectome.load(param, connName);

%% find soma for whole cells
somaAggloIds = find(conn.denMeta.targetClass == 'Somata');
somaAgglos = conn.dendrites(somaAggloIds);
somaLUT = Agglo.buildLUT(maxSegId, somaAgglos);

wcT = table;
wcT.id = reshape(1:numel(wc), [], 1);

% all segment IDs
wcT.allSegIds = arrayfun( ...
    @Agglo.fromSuperAgglo, ...
    wc, 'UniformOutput', false);

% axon segments
wcT.segIds = arrayfun( ...
    @(a) unique(a.nodes(a.axon ...
        & ~isnan(a.nodes(:, 4)), 4)), ...
    wc, 'UniformOutput', false);

% get rid of cells without axon segments
wcT(cellfun(@isempty, wcT.segIds), :) = [];

% find soma
wcT.somaId = cellfun( ...
    @(segIds) setdiff(somaLUT(segIds), 0), ...
    wcT.allSegIds, 'UniformOutput', false);
wcT.allSegIds = [];

wcT(~cellfun(@isscalar, wcT.somaId), :) = [];
wcT.somaId = cell2mat(wcT.somaId);

%% calculate distance to soma surface
% from `+connectEM/+Connectome/plotWholeCellInputs.m`
wcT.nodeDists = cell(size(wcT.segIds));

for curIdx = 1:size(wcT, 1)
    curAxonSegIds = wcT.segIds{curIdx};
    curSomaSegIds = somaAgglos{wcT.somaId(curIdx)};
    
    curSegIds = unique(cat(1, curAxonSegIds, curSomaSegIds));
   [~, curAxonSegIds] = ismember(curAxonSegIds, curSegIds);
   [~, curSomaSegIds] = ismember(curSomaSegIds, curSegIds);
    assert(all(curAxonSegIds) && all(curSomaSegIds));
    
    % build weights
    curDists = segCentroids(curSegIds, :);
    curDists = curDists .* param.raw.voxelSize;
    curDists = pdist(curDists);
    assert(all(curDists));
    
    curDists = squareform(curDists);
    
    % NOTE(amotta): Set distance between somatic segments to zero.
    %
    % True zeros cannot be used (not even with the additional 'Weights'
    % option of `graphminspantree`) since the output matrix will contain
    % zeros for both absent and somatic edges.
    curDists(curSomaSegIds, curSomaSegIds) = eps;
    
    curAdj = graphminspantree( ...
        sparse(curDists), curSomaSegIds(1));
    curAdj = logical(curAdj);
    
    % now we can use true zeros
    curDists(curSomaSegIds, curSomaSegIds) = 0;

    curDists = graphshortestpath( ...
        curAdj, curSomaSegIds(1), ...
        'Weights', curDists(curAdj), ...
        'Directed', false);
    
    curDists = curDists(curAxonSegIds);
    curDists = reshape(curDists, [], 1);
    wcT.nodeDists{curIdx} = curDists;
end

%% build presynapse table
preSynT = table;
preSynT.id = repelem( ...
    reshape(1:size(syn.synapses, 1), [], 1), ...
    cellfun(@numel, syn.synapses.presynId));
preSynT.segId = cell2mat(syn.synapses.presynId);

% assign dendrite IDs to synapses
dendLUT = Agglo.buildLUT(maxSegId, conn.dendrites);
syn.synapses.dendId = cellfun( ...
    @(ids) reshape(setdiff(dendLUT(ids), 0), [], 1), ...
    syn.synapses.postsynId, 'UniformOutput', false);

%% collect output synapses
wcT.synapses = cell(size(wcT.id));

for curIdx = 1:size(wcT, 1)
    curSegIds = wcT.segIds{curIdx};
    
    % find synapses
   [~, curSynIds] = ismember(curSegIds, preSynT.segId);
    curSynIds = setdiff(preSynT.id(setdiff(curSynIds, 0)), 0);
    
    % build table
    curSynT = table;
    curSynT.id = curSynIds;
    
    curSynT.segIdx = cellfun( ...
        @(ids) intersect(ids, wcT.segIds{curIdx}), ...
        syn.synapses.presynId(curSynT.id), 'UniformOutput', false);
    
    % assign synapse to largest segment
   [~, curSegIdx] = cellfun(@(ids) max(segMass(ids)), curSynT.segIdx);
    curSegIds = arrayfun(@(ids, idx) ids{1}(idx), curSynT.segIdx, curSegIdx);
    
    % translate to relative
   [~, curSynT.segIdx] = ismember(curSegIds, wcT.segIds{curIdx});
   
    % add dendrite IDs
    curSynReps = repelem( ...
        reshape(1:size(curSynT, 1), [], 1), ...
        cellfun(@numel, syn.synapses.dendId(curSynT.id)));
    
    curSynT = curSynT(curSynReps, :);
    curSynT.dendId = cell2mat(syn.synapses.dendId(curSynT.id));
    
    wcT.synapses{curIdx} = curSynT;
end

%% plotting
conn.denMeta.targetClass = reordercats(conn.denMeta.targetClass);
[targetClasses, ~, conn.denMeta.targetClassId] = unique(conn.denMeta.targetClass);

targetClasses(targetClasses == 'ApicalDendrite')     = 'AD';
targetClasses(targetClasses == 'SmoothDendrite')     = 'SD';
targetClasses(targetClasses == 'OtherDendrite')      = 'Other';
targetClasses(targetClasses == 'AxonInitialSegment') = 'AIS';
targetClasses(targetClasses == 'Somata')             = 'Soma';

targetClasses = arrayfun( ...
    @char, targetClasses, 'UniformOutput', false);

mkdir(outputDir);

for curIdx = 1:size(wcT, 1)
    curWcId = wcT.id(curIdx);
    curSyns = wcT.synapses{curIdx};
    if isempty(curSyns); continue; end
    
    curSyns.isSpine = syn.isSpineSyn(curSyns.id);
    curSyns.targetClassId = conn.denMeta.targetClassId(curSyns.dendId);
    
    curSyns.dist = wcT.nodeDists{curIdx}(curSyns.segIdx);
    curSyns.dist = curSyns.dist / 1E3;
    
    curMaxDist = 10 * ceil(max(curSyns.dist) / 10);
    curBinEdges = 0:10:curMaxDist;
    
    curPlot = @(ax, data) ...
        histogram( ...
            ax, data, ...
            'BinEdges', curBinEdges, ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2);
        
    curFig = figure('visible', 'off');
    curFig.Position(3:4) = [960, 1100];
    
    curAx = subplot(1 + numel(targetClasses), 1, 1);
    curAx.TickDir = 'out';
    hold(curAx, 'on');
    
	curPlot(curAx, curSyns.dist);
    curHist = curPlot(curAx, curSyns.dist(curSyns.isSpine));
    curHist.LineStyle = '--';
    
    xlim(curAx, curBinEdges([1, end]));
    ylabel(curAx, 'All');
    
    curAx.YLim(1) = 0;
    curYlim = curAx.YLim;
    
    legend(curAx, ...
        'All', 'Onto spines', ...
        'Location', 'West');
    
    title(curAx, { ...
        sprintf('Outputs from whole cell %d', curWcId); ...
        sprintf('%s (%s)', info.filename, info.git_repos{1}.hash)}, ...
        'FontWeight', 'normal', 'FontSize', 10);
    
    for curClassIdx = 1:numel(targetClasses)
        curAx = subplot( ...
            1 + numel(targetClasses), ...
            1, 1 + curClassIdx);
        curAx.TickDir = 'out';
        hold(curAx, 'on');
        
        curClassMask = (curSyns.targetClassId == curClassIdx);
        curPlot(curAx, curSyns.dist(curClassMask));
        
        curHist = curPlot(curAx, ...
            curSyns.dist(curClassMask & curSyns.isSpine));
        curHist.LineStyle = '--';
        
        ylim(curYlim);
        xlim(curAx, curBinEdges([1, end]));
        ylabel(curAx, targetClasses{curClassIdx});
    end
   
    xlabel(curAx, 'Distance to soma (µm)');
   
    curFigName = sprintf('output-distribution_whole-cell-%d.png', curWcId);
    export_fig('-r172', fullfile(outputDir, curFigName), curFig);
    clear curFig;
end