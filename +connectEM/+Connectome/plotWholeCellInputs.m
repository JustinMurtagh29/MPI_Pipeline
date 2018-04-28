% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

% Set output directory to write figures to disk instead of displaying them.
outputDir = '';

connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');
wcFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto-and-manual.mat');
somaFile  = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_03_v2.mat');

% debugging
debugDir = '';
debugCellIds = [];

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

segMass = Seg.Global.getSegToSizeMap(param);

[conn, syn] = connectEM.Connectome.load(param, connFile);

% Complete dendrite meta information
denMeta = load(conn.info.param.dendriteFile);
denMetaWcIdx = denMeta.idxWholeCells(conn.denMeta.parentId);
denMetaSomaIdx = denMeta.idxSomata(conn.denMeta.parentId);

assert(~any(denMetaWcIdx & denMetaSomaIdx));
conn.denMeta.wcId = max(denMetaWcIdx, denMetaSomaIdx);

wcData = load(wcFile);
somaData = load(somaFile);

%% split axons into exc. and inh.
%{
conn.axonMeta.spineSynFrac = ...
    conn.axonMeta.spineSynCount ...
 ./ conn.axonMeta.synCount;

conn.axonMeta.isExc = (conn.axonMeta.priSpineSynFrac >= 0.5);
conn.axonMeta.isInh = ~conn.axonMeta.isExc;
%}

%% Complete whole cell somata
wcT = table;
wcT.id = wcData.idxWholeCells(wcData.indWholeCells);
wcT.agglo = wcData.dendrites(wcData.indWholeCells);
SuperAgglo.check(wcT.agglo);

% Find corresponding somata
[~, somaIds] = ismember(wcT.id, somaData.idxSomata);
wcT.somaAgglo = somaData.dendrites(somaIds);

% NOTE(amotta): There is a bug in the soma super-agglomerates which allows
% them to be disconnected. Let's fix this by introducing random edges. This
% is not a problem since we do not care about distances within the soma.
wcT.somaAgglo = SuperAgglo.connect(wcT.somaAgglo);
SuperAgglo.check(wcT.somaAgglo);

% Merge somata with whole cells
wcT.agglo = arrayfun(@SuperAgglo.merge, wcT.agglo, wcT.somaAgglo);
wcT.agglo = SuperAgglo.clean(wcT.agglo);

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

%% collect input synapses
wcT.synapses = cell(size(wcT.id));

for curIdx = 1:size(wcT, 1)
    curAgglo = wcT.agglo(curIdx);
    
    curPostIds = conn.denMeta.wcId == wcT.id(curIdx);
    curPostIds = conn.denMeta.id(curPostIds);
    
    curConnMask = ismember( ...
        conn.connectome.edges(:, 2), curPostIds);
    curConnRows = conn.connectome(curConnMask, :);
    curAreaRows = conn.connectomeMeta.contactArea(curConnMask);
    
    curSynT = table;
    curSynT.id = cell2mat(curConnRows.synIdx);
    curSynT.area = cell2mat(curAreaRows);
    
    curSynT.nodeId = cellfun( ...
        @(ids) intersect(ids, curAgglo.nodes(:, 4)), ...
        syn.synapses.postsynId(curSynT.id), ...
        'UniformOutput', false);
    
    % assign synapse to largest segment
   [~, curNodeId] = cellfun( ...
        @(ids) max(segMass(ids)), curSynT.nodeId);
    curNodeIds = arrayfun( ...
        @(ids, idx) ids{1}(idx), curSynT.nodeId, curNodeId);
    
    % translate to relative 
   [~, curSynT.curNodeId] = ismember( ...
       double(curNodeIds), curAgglo.nodes(:, 4));
    
    curSynT.axonId = repelem( ...
        curConnRows.edges(:, 1), ...
        cellfun(@numel, curConnRows.synIdx));
    
    wcT.synapses{curIdx} = curSynT;
end

%% debugging
if ~isempty(debugDir)
    skel = skeleton();
    skel = Skeleton.setParams4Pipeline(skel, param);
    skel = skel.setDescription(sprintf( ...
        '%s (%s)', info.filename, info.git_repos{1}.hash));

    for curIdx = 1:numel(debugCellIds)
        curId = debugCellIds(curIdx);
        
        curEdges = wcT.edges{curId};
        curNodeIds = wcT.segIds{curId};
        curDistsUm = wcT.nodeDists{curId} / 1E3;
        curSynIdx = wcT.synapses{curId}.segIdx;

        % generate skeleton
        curSkel = skel.addTree( ...
            sprintf('Whole cell %d', curId), ...
            ceil(segCentroids(curNodeIds, :)), curEdges);

        % add distance annotations
        curComments = arrayfun( ...
            @(segId, distUm) sprintf( ...
                'Synapse. Segment %d. %.1f µm to soma.', segId, distUm), ...
            curNodeIds(curSynIdx), curDistsUm(curSynIdx), 'UniformOutput', false);
       [curSkel.nodesAsStruct{end}(curSynIdx).comment] = deal(curComments{:});

        % add branchpoints
        curSkel = curSkel.addBranchpoint(curSynIdx);

        curSkelName = sprintf( ...
            '%0*d_whole-cell-%d.nml', ...
            ceil(log10(1 + size(wcT, 1))), curIdx, curId);
        curSkel.write(fullfile(debugDir, curSkelName));
    end
end

%% plotting
for curIdx = 1:size(wcT, 1)
    curSyns = wcT.synapses{curIdx};
    if isempty(curSyns); continue; end
    
    curFig = figure('visible', 'off');
    curFig.Position(3:4) = [1075, 600];
    
    curSyns.isSpine = syn.isSpineSyn(curSyns.id);
    curSyns.isSoma = logical(somaLUT( ...
        wcT.segIds{curIdx}(curSyns.segIdx)));
    
    curSyns.isExc = conn.axonMeta.isExc(curSyns.axonId);
    curSyns.isInh = conn.axonMeta.isInh(curSyns.axonId);
    curSyns.isTc = conn.axonMeta.isThalamocortical(curSyns.axonId);
    
    curSyns.dist = wcT.nodeDists{curIdx}(curSyns.segIdx);
    curSyns.dist = curSyns.dist / 1E3;
    
    % move soma synapses to separate bin
    curSyns.dist(curSyns.isSoma) = -eps;
    
    % find repeated synapses
    curRepT = table;
   [curRepT.axonId, ~, curSynCount] = unique(curSyns.axonId);
    curRepT.synDist = accumarray( ...
        curSynCount, curSyns.dist, [], ...
        @(dist) {sort(dist, 'ascend')});
    
    curRepT.synCount = accumarray(curSynCount, 1);
    curRepT(curRepT.synCount < 2, :) = [];
    clear curSynCount;
    
    % sorting
   [~, curSortIds] = sort(cellfun( ...
       @median, curRepT.synDist), 'descend');
    curRepT = curRepT(curSortIds, :);
    
    curMaxDist = 10 * ceil(max(curSyns.dist) / 10);
    curBinEdges = -10:10:curMaxDist;
    
    curPlot = @(ax, data) ...
        histogram( ...
            ax, data, ...
            'BinEdges', curBinEdges, ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2);
        
    curAx = subplot(4, 1, 1);
    hold(curAx, 'on');
    
    curPlot(curAx, curSyns.dist);
    curPlot(curAx, curSyns.dist(curSyns.isSpine));
    curPlot(curAx, curSyns.dist(curSyns.isSoma));
    
    curAx.TickDir = 'out';
    curAx.Position(3) = 0.8 - curAx.Position(1);
    
    xlim(curAx, curBinEdges([1, end]));
    ylabel(curAx, 'Synapses');
    
    title(curAx, { ...
        sprintf('Inputs onto whole cell %d', curIdx); ...
        sprintf('%s (%s)', info.filename, info.git_repos{1}.hash)}, ...
        'FontWeight', 'normal', 'FontSize', 10);
    
    curLeg = legend(curAx, ...
        'All', 'Onto spines', 'Onto soma', ...
        'Location', 'EastOutside');
    curLeg.Position([1, 3]) = [0.82, (0.98 - 0.82)];
    
    curAx = subplot(4, 1, 2);
    hold(curAx, 'on');
    
    curDiscretize = ...
        @(data) accumarray( ...
            discretize(data, curBinEdges), ...
            1, [numel(curBinEdges) - 1, 1]);
    curBinsAll = curDiscretize(curSyns.dist);
    
    curPlot = ...
        @(ax, data) histogram(ax, ...
            'BinCount', min( ...
                curDiscretize(data) ...
                ./ curBinsAll, 1), ...
            'BinEdges', curBinEdges, ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2);
        
    curPlot(curAx, curSyns.dist(curSyns.isExc));
    curPlot(curAx, curSyns.dist(curSyns.isTc));
    
    if ~isempty(curRepT)
        % `histogram` fails if there are no synapses
        curPlot(curAx, cell2mat(curRepT.synDist));
    end
    
    curAx.TickDir = 'out';
    curAx.Position(3) = 0.8 - curAx.Position(1);
    
    xlim(curAx, curBinEdges([1, end]));
    ylabel(curAx, 'Fraction of synapses');
    ylim(curAx, [0, 1]);
    
    curLeg = legend(curAx, ...
        'Excitatory', ...
        'Thalamocortical', ...
        'Repeated', ...
        'Location', 'EastOutside');
    curLeg.Position([1, 3]) = [0.82, (0.98 - 0.82)];
    
    % plot synaptic clustering
    curAx = subplot(4, 1, 3);
    hold(curAx, 'on');
    
    for curRepIdx = 1:size(curRepT, 1)
        % randomize soma synapse locations
        curRepSynDists = curRepT.synDist{curRepIdx};
        curRepSynDists(curRepSynDists < 0) = curBinEdges(1) + ...
            diff(curBinEdges(1:2)) * rand(sum(curRepSynDists < 0), 1);
        
        plot( ...
            curAx, curRepSynDists, ...
            repelem(curRepIdx, curRepT.synCount(curRepIdx)), ...
            'k.-', 'MarkerSize', 12);
    end
    
    curAx.TickDir = 'out';
    curAx.Position(3) = 0.8 - curAx.Position(1);
    
    xlim(curAx, curBinEdges([1, end]));
    ylim(curAx, [0, max(1, size(curRepT, 1))]);
    yticks(curAx, curAx.YLim);
    ylabel(curAx, 'Axons');
    
    % plot (spine) synapse sizes
    curAx = subplot(4, 1, 4);
    
    curSpineSyns = curSyns(curSyns.isSpine, :);
    curSpineSyns.binId = discretize( ...
        curSpineSyns.dist, curBinEdges);
    curSpineSynArea = accumarray( ...
        curSpineSyns.binId, curSpineSyns.area, ...
       [numel(curBinEdges) - 1, 1], @median, 0);
   
    histogram(curAx, ...
        'BinCount', curSpineSynArea, ...
        'BinEdges', curBinEdges, ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2);
    
    curAx.TickDir = 'out';
    curAx.Position(3) = 0.8 - curAx.Position(1);
    
    ylabel(curAx, 'Median ASI (µm²)');
    xlabel(curAx, 'Distance to soma (µm)');
    xlim(curAx, curBinEdges([1, end]));
    
    % save figure
    if ~isempty(outputDir)
        curFigFile = fullfile(outputDir, sprintf( ...
            'input-distribution_whole-cell-%d.png', curIdx));
        export_fig(curFigFile, curFig);
        clear curFig;
    end
end

%% load cached results
%{
% TODO(amotta): remove this
wcT = load(fullfile(outputDir, 'whole-cell-table.mat'));
wcT = wcT.wcT;
%}

%% lump inputs
% build synapse table and fix `segIdx`
segIdxOff = cumsum(cellfun(@numel, wcT.segIds));
segIdxOff = [0; segIdxOff(1:(end - 1))];

lumpedSynapses = cat(1, wcT.synapses{:});
lumpedSynapses.segIdx = lumpedSynapses.segIdx ...
    + repelem(segIdxOff, cellfun(@height, wcT.synapses));

lumpedWcT = struct;
lumpedWcT.segIds = cell2mat(wcT.segIds);
lumpedWcT.nodeDists = cell2mat(wcT.nodeDists);
lumpedWcT.synapses = lumpedSynapses;

%% plot lumped inputs
curFig = figure();
curFig.Position(3:4) = [860, 480];

curSyns = lumpedWcT.synapses;
curSyns.isSpine = syn.isSpineSyn(curSyns.id);
curSyns.isSoma = logical(somaLUT(lumpedWcT.segIds(curSyns.segIdx)));

curSyns.isExc = conn.axonMeta.isExc(curSyns.axonId);
curSyns.isInh = conn.axonMeta.isInh(curSyns.axonId);
curSyns.isTc  = conn.axonMeta.isThalamocortical(curSyns.axonId);

curSyns.dist = lumpedWcT.nodeDists(curSyns.segIdx);
curSyns.dist = curSyns.dist / 1E3;

% move soma synapses to separate bin
curSyns.dist(curSyns.isSoma) = -eps;

curMaxDist = 10 * ceil(max(curSyns.dist) / 10);
curBinEdges = -10:10:curMaxDist;

curPlot = @(ax, data) ...
    histogram( ...
        ax, data, ...
        'BinEdges', curBinEdges, ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2);

curAx = subplot(3, 1, 1);
hold(curAx, 'on');

curPlot(curAx, curSyns.dist);
curPlot(curAx, curSyns.dist(curSyns.isSpine));
curPlot(curAx, curSyns.dist(curSyns.isSoma));

curAx.TickDir = 'out';
curAx.Position(3) = 0.8 - curAx.Position(1);

xlim(curAx, curBinEdges([1, end]));
ylabel(curAx, 'Synapses');

title(curAx, { ...
    'Lumped inputs onto soma-based reconstructions'; ...
    sprintf('%s (%s)', info.filename, info.git_repos{1}.hash)}, ...
    'FontWeight', 'normal', 'FontSize', 10);

curLeg = legend(curAx, ...
    'All', 'Onto spines', 'Onto soma', ...
    'Location', 'EastOutside');
curLeg.Position([1, 3]) = [0.82, (0.98 - 0.82)];

% plot fraction of axon types
curAx = subplot(3, 1, 2);
hold(curAx, 'on');

curDiscretize = ...
    @(data) accumarray( ...
        discretize(data, curBinEdges), ...
        1, [numel(curBinEdges) - 1, 1]);
curBinAll = curDiscretize(curSyns.dist);

curSpineFrac = curDiscretize(curSyns.dist(curSyns.isSpine));
curExcFrac = curDiscretize(curSyns.dist(curSyns.isExc));
curTcFrac = curDiscretize(curSyns.dist(curSyns.isTc));

curPlot = ...
    @(ax, data) histogram(ax, ...
        'BinEdges', curBinEdges, ...
        'BinCounts', data, ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2);

curPlot(curAx, curExcFrac ./ curBinAll);
curPlot(curAx, curSpineFrac ./ curBinAll);
curPlot(curAx, curTcFrac ./ curBinAll);

curLeg = legend(curAx, ...
    'Excitatory', ...
    'Spine', ...
    'Thalamocortical', ...
    'Location', 'EastOutside');
curLeg.Position([1, 3]) = [0.82, (0.98 - 0.82)];

curAx.TickDir = 'out';
curAx.Position(3) = 0.8 - curAx.Position(1);

xlim(curAx, curBinEdges([1, end]));
ylabel(curAx, 'Fraction of synapses');

% plot axon types
curAx = subplot(3, 1, 3);
hold(curAx, 'on');

curPlot = @(ax, data) ...
    histogram( ...
        ax, data, ...
        'Normalization', 'probability', ...
        'BinEdges', curBinEdges, ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2);

curPlot(curAx, curSyns.dist(curSyns.isExc));
curPlot(curAx, curSyns.dist(curSyns.isInh));
curPlot(curAx, curSyns.dist(curSyns.isTc));

curAx.TickDir = 'out';
curAx.Position(3) = 0.8 - curAx.Position(1);

xlim(curAx, curBinEdges([1, end]));
ylabel(curAx, 'Synapse probability');

curLeg = legend(curAx, ...
    'Excitatory', ...
    'Inhibitory', ...
    'Thalamocortical', ...
    'Location', 'EastOutside');
curLeg.Position([1, 3]) = [0.82, (0.98 - 0.82)];

xlabel(curAx, 'Distance to soma (µm)');

%% calculate fraction of spine synapses from TC axons
allWcSyns = cat(1, wcT.synapses{:});

% only consider spine synapses
allWcSyns.isSpine = syn.isSpineSyn(allWcSyns.id);
allWcSyns(~allWcSyns.isSpine, :) = [];

% check if from TC axon
allWcSyns.isTc = conn.axonMeta.isThalamocortical(allWcSyns.axonId);
allWcTcSpineFrac = mean(allWcSyns.isTc);

fprintf( ...
   ['Fraction of spine synapses ', ...
    'from TC axons: %.1f %%\n'], ...
    100 * allWcTcSpineFrac);
