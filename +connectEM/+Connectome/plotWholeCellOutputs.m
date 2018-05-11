% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
outputDir = '/home/amotta/Desktop/whole-cell-output-distributions';

wcFile = fullfile(rootDir, 'aggloState', 'wholeCells_GTAxon_08_v4.mat');
somaFile  = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_03_v2.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

[conn, syn] = connectEM.Connectome.load(param, connFile);

wcData = load(wcFile);
somaData = load(somaFile);

maxSegId = Seg.Global.getMaxSegId(param);
segSizes = Seg.Global.getSegToSizeMap(param);

%% Build whole cell table
wcAgglos = SuperAgglo.clean(wcData.wholeCells);
[wcAgglos.axon] = deal(wcData.wholeCells.axon);
wcSegEqs = Agglo.fromSuperAgglo(wcAgglos);

wcCellIds = ...
    conn.denMeta.targetClass == 'WholeCell';
wcCellIds = Agglo.buildLUT( ...
    maxSegId, ...
    conn.dendrites(wcCellIds), ...
    conn.denMeta.cellId(wcCellIds));
wcCellIds = cellfun(@(segIds) mode( ...
    nonzeros(wcCellIds(segIds))), wcSegEqs);

wcT = table;
wcT.id = wcCellIds;
wcT.agglo = wcAgglos;

% Fix and prepare soma super-agglomerates
somaAgglos = somaData.dendrites(somaData.indSomata);
somaAgglos = SuperAgglo.connect(somaAgglos);
SuperAgglo.check(somaAgglos);

% Find soma agglomerate
somaLUT = Agglo.buildLUT( ...
    maxSegId, Agglo.fromSuperAgglo(somaAgglos));
somaIds = cellfun(@(segIds) mode( ...
    nonzeros(somaLUT(segIds))), wcSegEqs);
assert(numel(somaIds) == numel(unique(somaIds)));

% Fix soma agglomerate
wcT.somaAgglo = somaAgglos(somaIds);
wcT.agglo = arrayfun(@SuperAgglo.merge, wcT.agglo, wcT.somaAgglo);

% Axon nodes by size
wcT.axonNodeIds = arrayfun( ...
    @(wc) reshape(find(wc.axon), [], 1), ...
    wcAgglos, 'UniformOutput', false);

% Add meta data for plotting
wcT.title = arrayfun( ...
    @(id) sprintf('whole cell %d', id), ...
    wcT.id, 'UniformOutput', false);
wcT.tag = strrep(wcT.title, ' ', '-');

%% Calculate distance to soma surface
wcT.nodeDists = cell(size(wcT.id));

for curIdx = 1:size(wcT, 1)
    curAgglo = wcT.agglo(curIdx);
    curNodeCount = size(curAgglo.nodes, 1);
    
    curSomaSegIds = Agglo.fromSuperAgglo(wcT.somaAgglo(curIdx));
    curSomaNodeMask = ismember(curAgglo.nodes(:, 4), curSomaSegIds);
    
    curSomaNodeId = find(curSomaNodeMask, 1);
    assert(~isempty(curSomaNodeId));
    
    curDists = ...
        curAgglo.nodes(curAgglo.edges(:, 1), 1:3) ...
      - curAgglo.nodes(curAgglo.edges(:, 2), 1:3);
    curDists = curDists .* param.raw.voxelSize;
    curDists = sqrt(sum(curDists .* curDists, 2));
    assert(all(curDists));
    
    % Travelling within soma is for free!
    curDists(all(curSomaNodeMask(curAgglo.edges), 2)) = 0;
    
    curAdj = sparse( ...
        curAgglo.edges(:, 2), curAgglo.edges(:, 1), ...
        true, curNodeCount, curNodeCount);

    curDists = graphshortestpath( ...
        curAdj, curSomaNodeId, ...
        'Weights', curDists, ...
        'Directed', false);
    
    % Sanity checks
    assert(~any(curDists(curSomaNodeMask)));
    assert(all(curDists(~curSomaNodeMask)));
    
    curDists = reshape(curDists, [], 1);
    wcT.nodeDists{curIdx} = curDists;
end

%% Build presynapse table
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

%% Collect output synapses
last = @(vals) vals(end);
wcT.synapses = cell(size(wcT.id));

for curIdx = 1:size(wcT, 1)
    curCellId = wcT.id(curIdx);
    curSegIds = wcT.agglo(curIdx).nodes(:, 4);
    
    curAxonNodeIds = wcT.axonNodeIds{curIdx};
    curAxonNodeIds(isnan(curSegIds(curAxonNodeIds))) = [];
    
    % NOTE(amotta): Sort axon nodes by increasing size. This will make it
    % much easier to assign synapses to the largest presynaptic segment.
    curAxonNodeIds = Util.sortBy( ...
        curAxonNodeIds, segSizes(curSegIds(curAxonNodeIds)));
    curAxonSegIds = curSegIds(curAxonNodeIds);
    
    % Find synapses
   [~, curSynIds] = ismember(curAxonSegIds, preSynT.segId);
    curSynIds = unique(preSynT.id(setdiff(curSynIds, 0)));
    
    % Build table
    curSynT = table;
    curSynT.id = curSynIds;
    
   [~, curSynT.nodeId] = cellfun( ...
        @(ids) intersect(curAxonSegIds, ids, 'stable'), ...
        syn.synapses.presynId(curSynT.id), 'UniformOutput', false);
    curSynT.nodeId = cellfun( ...
        @(ids) curAxonNodeIds(ids(end)), curSynT.nodeId);
    
    % Add dendrite IDs
    curSynReps = repelem( ...
        reshape(1:size(curSynT, 1), [], 1), ...
        cellfun(@numel, syn.synapses.dendId(curSynT.id)));
    
    curSynT = curSynT(curSynReps, :);
    curSynT.dendId = cell2mat(syn.synapses.dendId(curSynT.id));
    
    % Discard false positive synapses from the axon onto soma and / or axon
    % initial segment of the same cell.
    curOwnDendIds = conn.denMeta.id(conn.denMeta.cellId == curCellId);
    curSynT(ismember(curSynT.dendId, curOwnDendIds), :) = [];
    
    wcT.synapses{curIdx} = curSynT;
end

%% Plot average over all axons with 10 or more synapses
wcGroups = struct;
wcGroups(1).wcIds = find( ...
    cellfun(@height, wcT.synapses) >= 10);
wcGroups(1).title = sprintf( ...
    'axons with ≥ 10 synapses (n = %d)', ...
    numel(wcGroups(1).wcIds));
wcGroups(1).tag = 'min-10-syn';

%% Generate extended whole cell table with queen axons
extWcT = wcT;
for curIdx = 1:numel(wcGroups)
    curWcGroup = wcGroups(curIdx);
    curWcT = extWcT(curWcGroup.wcIds, :);
    
    nodeIdOff = cumsum(cellfun(@numel, curWcT.nodeDists));
    nodeIdOff = [0; nodeIdOff(1:(end - 1))];
    
    axonNodeIds = cat(1, curWcT.axonNodeIds{:});
    axonNodeIds = axonNodeIds + repelem(nodeIdOff, ...
        cellfun(@(ids) size(ids, 1), curWcT.axonNodeIds));

    lumpedSynapses = cat(1, curWcT.synapses{:});
    lumpedSynapses.nodeId = lumpedSynapses.nodeId ...
        + repelem(nodeIdOff, cellfun(@height, curWcT.synapses));

    curQueenWcT = table;
    curQueenWcT.id = 0;
    curQueenWcT.axonNodeIds = {axonNodeIds};

    curQueenWcT.agglo = struct;
    curQueenWcT.agglo.nodes = cat(1, curWcT.agglo.nodes);
    curQueenWcT.agglo.edges = zeros(0, 2);

    curQueenWcT.somaAgglo = struct;
    curQueenWcT.somaAgglo.nodes = cat(1, curWcT.somaAgglo.nodes);
    curQueenWcT.somaAgglo.edges = zeros(0, 2);

    curQueenWcT.nodeDists = {cell2mat(curWcT.nodeDists)};
    curQueenWcT.synapses = {lumpedSynapses};
    
    curQueenWcT.title = curWcGroup.title;
    curQueenWcT.tag = curWcGroup.tag;

    extWcT = cat(1, extWcT, curQueenWcT);
end

%% Determine bin edges
binSize = 20;

% Find farthest synapse
binEdges = cellfun( ...
    @(syns, nodeDists) max(nodeDists(syns.nodeId)), ...
    extWcT.synapses, extWcT.nodeDists, 'UniformOutput', false);
binEdges(cellfun(@isempty, binEdges)) = [];
binEdges = max(cell2mat(binEdges)) / 1E3;

% Build bin edges
binEdges = binSize * ceil(binEdges / binSize);
binEdges = 0:binSize:binEdges;

%% Plotting
[targetClasses, ~, conn.denMeta.targetClassId] = unique(conn.denMeta.targetClass);

targetClasses(targetClasses == 'ApicalDendrite')     = 'AD';
targetClasses(targetClasses == 'SmoothDendrite')     = 'SD';
targetClasses(targetClasses == 'OtherDendrite')      = 'Other';
targetClasses(targetClasses == 'AxonInitialSegment') = 'AIS';
targetClasses(targetClasses == 'Somata')             = 'Soma';

targetClasses = arrayfun( ...
    @char, targetClasses, 'UniformOutput', false);

for curIdx = 1:size(extWcT, 1)
    curTag = extWcT.tag{curIdx};
    curTitle = extWcT.title{curIdx};
    curSyns = extWcT.synapses{curIdx};
    if isempty(curSyns); continue; end
    
    curSyns.isSpine = syn.isSpineSyn(curSyns.id);
    curSyns.targetClassId = conn.denMeta.targetClassId(curSyns.dendId);
    
    curSyns.dist = extWcT.nodeDists{curIdx}(curSyns.nodeId);
    curSyns.dist = curSyns.dist / 1E3;
    
    curPlot = @(ax, data) ...
        histogram( ...
            ax, data, ...
            'BinEdges', binEdges, ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2, ...
            'FaceAlpha', 1);
    
    if ~isempty(outputDir)
        curFig = figure('visible', 'off');
    else
        curFig = figure();
    end
    
    curFig.Color = 'white';
    curFig.Position(3:4) = [960, 1100];
    
    curAx = subplot(1 + numel(targetClasses), 1, 1);
    curAx.TickDir = 'out';
    hold(curAx, 'on');
    
    curHist = curPlot(curAx, curSyns.dist(curSyns.isSpine));
    curHist.EdgeColor = 'magenta';
    
    curHist = curPlot(curAx, curSyns.dist(~curSyns.isSpine));
    curHist.EdgeColor = 'black';
    
    xlim(curAx, binEdges([1, end]));
    ylabel(curAx, 'All');
    
    curAx.YLim(1) = 0;
    curYlim = curAx.YLim;
    
    title(curAx, { ...
        sprintf('Outputs from %s', curTitle); ...
        sprintf('%s (%s)', info.filename, info.git_repos{1}.hash)}, ...
        'FontWeight', 'normal', 'FontSize', 10);
    
    for curClassIdx = 1:numel(targetClasses)
        curAx = subplot( ...
            1 + numel(targetClasses), ...
            1, 1 + curClassIdx);
        curAx.TickDir = 'out';
        hold(curAx, 'on');
        
        curClassMask = (curSyns.targetClassId == curClassIdx);
        
        curHist = curPlot(curAx, ...
            curSyns.dist(curClassMask &  curSyns.isSpine));
        curHist.EdgeColor = 'magenta';
        
        curHist = curPlot(curAx, ...
            curSyns.dist(curClassMask & ~curSyns.isSpine));
        curHist.EdgeColor = 'black';
        
        ylim(curYlim);
        xlim(curAx, binEdges([1, end]));
        ylabel(curAx, targetClasses{curClassIdx});
    end
   
    xlabel(curAx, 'Distance to soma (µm)');
   
    if ~isempty(outputDir)
        curFigFile = sprintf('output-distribution_%s', curTag);
        curFigFile = fullfile(outputDir, curFigFile);
        
        export_fig(strcat(curFigFile, '.png'), curFig);
        export_fig(strcat(curFigFile, '.eps'), curFig);
        clear curFig;
    end
end
