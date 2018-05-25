% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

% Set output directory to write figures to disk instead of displaying them.
plotDir = '';
plotShow = false;

connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');
wcFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto-and-manual.mat');
somaFile  = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_03_v2.mat');

splitNmlDir = fileparts(fileparts(mfilename('fullpath')));
splitNmlDir = fullfile(splitNmlDir, '+WholeCell', '+Script', 'annotations');

% debugging
debugDir = '';
debugCellIds = [];

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

segMass = Seg.Global.getSegToSizeMap(param);
segPoint = Seg.Global.getSegToPointMap(param);

[conn, syn] = connectEM.Connectome.load(param, connFile);

wcData = load(wcFile);
somaData = load(somaFile);

%% NML files for whole cell splitting
splitNmlT = connectEM.WholeCell.loadSplitNmls(splitNmlDir);
splitNmlT.cellId = wcData.idxWholeCells(splitNmlT.aggloId);
assert(all(splitNmlT.cellId));

%% Split axons into exc. and inh.
conn.axonMeta.fullPriSpineSynFrac = ...
    conn.axonMeta.fullPriSpineSynCount ...
 ./ conn.axonMeta.fullSynCount;

conn.axonMeta.isExc = (conn.axonMeta.fullPriSpineSynFrac >= 0.5);
conn.axonMeta.isInh = ~conn.axonMeta.isExc;

%% Complete whole cell somata
wcT = table;
wcT.id = wcData.idxWholeCells(wcData.indWholeCells);
wcT.agglo = wcData.dendrites(wcData.indWholeCells);
SuperAgglo.check(wcT.agglo);

% Find corresponding somata
[~, somaIds] = ismember(wcT.id, somaData.idxSomata);
wcT.somaAgglo = somaData.dendrites(somaIds);

calcSomaPos = @(n) ...
    sum(segMass(n(:, 4)) .* n(:, 1:3), 1) ....
    ./ sum(segMass(n(:, 4)));
wcT.somaPos = cell2mat(arrayfun( ...
    @(a) calcSomaPos(a.nodes), ...
    wcT.somaAgglo, 'UniformOutput', false));

% NOTE(amotta): There is a bug in the soma super-agglomerates which allows
% them to be disconnected. Let's fix this by introducing random edges. This
% is not a problem since we do not care about distances within the soma.
wcT.somaAgglo = SuperAgglo.connect(wcT.somaAgglo);
SuperAgglo.check(wcT.somaAgglo);

% Merge somata with whole cells
wcT.agglo = arrayfun(@SuperAgglo.merge, wcT.agglo, wcT.somaAgglo);
wcT.agglo = SuperAgglo.clean(wcT.agglo);

wcT.title = arrayfun( ...
    @(id) sprintf('whole cell %d', id), ...
    wcT.id, 'UniformOutput', false);
wcT.tag = strrep(wcT.title, ' ', '-');

%% Calculate node distances
%  * to soma surface, and
%  * orthogonal to soma
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
    
    % Orthogonal distance
    curOrthoDists = curAgglo.nodes(:, 1:3) - wcT.somaPos(curIdx, :);
    curOrthoDists = curOrthoDists .* param.raw.voxelSize;
    
    wcT.nodeDists{curIdx} = curDists;
    wcT.nodeOrthoDists{curIdx} = curOrthoDists;
end

%% Collect input synapses
wcT.synapses = cell(size(wcT.id));

for curIdx = 1:size(wcT, 1)
    curAgglo = wcT.agglo(curIdx);
    
    curPostIds = ...
        conn.denMeta.cellId == wcT.id(curIdx) ...
      & ismember(conn.denMeta.targetClass, {'Somata', 'WholeCell'});
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
   [~, curSynT.nodeId] = ismember( ...
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

%% Try to find border cells
somaPos = cell2mat(arrayfun( ...
    @(s) mean(s.nodes(:, 1:3), 1), ...
    wcT.somaAgglo, 'UniformOutput', false));

voxelSize = param.raw.voxelSize;
somaDistXY = voxelSize(3) * min(pdist2( ...
    somaPos(:, 3), transpose(param.bbox(3, :))), [], 2);
somaDistXZ = voxelSize(2) * min(pdist2( ...
    somaPos(:, 2), transpose(param.bbox(2, :))), [], 2);
somaDistYZ = voxelSize(1) * min(pdist2( ...
    somaPos(:, 1), transpose(param.bbox(1, :))), [], 2);

yzWcIds = find(somaDistYZ >= 10E3 & ...
    (somaDistXY < 10E3 | somaDistXZ < 10E3));

wcGroups = struct;
wcGroups(1).wcIds = yzWcIds;
wcGroups(1).title = sprintf( ...
    'whole cells in YZ view (n = %d)', ...
    numel(wcGroups(1).wcIds));
wcGroups(1).tag = 'yz-view';

wcGroups(2).wcIds = wcT.id( ...
    cellfun(@height, wcT.synapses) > 100);
wcGroups(2).title = sprintf( ...
    'whole cells with > 100 synapses (n = %d)', ...
    numel(wcGroups(2).wcIds));
wcGroups(2).tag = '100-syn';

wcGroups(3).wcIds = wcT.id( ...
    cellfun(@height, wcT.synapses) > 1);
wcGroups(3).title = sprintf( ...
    'whole cells with > 1 synapses (n = %d)', ...
    numel(wcGroups(3).wcIds));
wcGroups(3).tag = '1-syn';

%% Generate queen neuron
extWcT = wcT;
for curIdx = 1:numel(wcGroups)
    curWcGroup = wcGroups(curIdx);
    curWcT = extWcT(curWcGroup.wcIds, :);
    
    nodeIdOff = cumsum(cellfun(@numel, curWcT.nodeDists));
    nodeIdOff = [0; nodeIdOff(1:(end - 1))];

    lumpedSynapses = cat(1, curWcT.synapses{:});
    lumpedSynapses.nodeId = lumpedSynapses.nodeId ...
        + repelem(nodeIdOff, cellfun(@height, curWcT.synapses));

    curQueenWcT = table;
    curQueenWcT.id = 0;

    curQueenWcT.agglo = struct;
    curQueenWcT.agglo.nodes = cat(1, curWcT.agglo.nodes);
    curQueenWcT.agglo.edges = zeros(0, 2);

    curQueenWcT.somaAgglo = struct;
    curQueenWcT.somaAgglo.nodes = cat(1, curWcT.somaAgglo.nodes);
    curQueenWcT.somaAgglo.edges = zeros(0, 2);
    
    curQueenWcT.somaPos = nan(1, 3);

    curQueenWcT.nodeDists = {cell2mat(curWcT.nodeDists)};
    curQueenWcT.nodeOrthoDists = {cell2mat(curWcT.nodeOrthoDists)};
    curQueenWcT.synapses = {lumpedSynapses};
    
    curQueenWcT.title = curWcGroup.title;
    curQueenWcT.tag = curWcGroup.tag;

    extWcT = cat(1, extWcT, curQueenWcT);
end

%% Plotting
for curIdx = size(extWcT, 1)
    curSyns = extWcT.synapses{curIdx};
    if isempty(curSyns); continue; end
    
    if plotShow
        curFig = figure(); %#ok
    elseif ~isempty(plotDir)
        curFig = figure('visible', 'off');
    else
        % No need to build plot
        continue;
    end
    
    curFig.Color = 'white';
    curFig.Position(3:4) = [1120, 660];
    
    curSyns.isSpine = syn.isSpineSyn(curSyns.id);
    curSyns.isSoma = ismember( ...
        extWcT.agglo(curIdx).nodes(curSyns.nodeId, 4), ...
        extWcT.somaAgglo(curIdx).nodes(:, 4));
    
    curSyns.isExc = conn.axonMeta.isExc(curSyns.axonId);
    curSyns.isInh = conn.axonMeta.isInh(curSyns.axonId);
    
    curSyns.axonClass = conn.axonMeta.axonClass(curSyns.axonId);
    
    curSyns.dist = extWcT.nodeDists{curIdx}(curSyns.nodeId);
    curSyns.dist = curSyns.dist / 1E3;
    
    % move soma synapses to separate bin
    curSyns.dist(curSyns.isSoma) = -eps;
    
    curMaxDist = 10 * ceil(max(curSyns.dist) / 10);
    curBinEdges = -5:5:curMaxDist;
    
    curPlot = @(ax, data) ...
        histogram( ...
            ax, data, ...
            'BinEdges', curBinEdges, ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2, ...
            'FaceAlpha', 1);
        
    curAx = subplot(4, 1, 1);
    hold(curAx, 'on');
    
    curPlot(curAx, curSyns.dist);
    curPlot(curAx, curSyns.dist(curSyns.isSpine));
    curPlot(curAx, curSyns.dist(curSyns.isSoma | ~curSyns.isSpine));
    
    curAx.TickDir = 'out';
    curAx.Position(3) = 0.8 - curAx.Position(1);
    
    xlim(curAx, curBinEdges([1, end]));
    ylabel(curAx, 'Synapses');
    
    title(curAx, { ...
        sprintf('Inputs onto %s', extWcT.title{curIdx}); ...
        sprintf('%s (%s)', info.filename, info.git_repos{1}.hash)}, ...
        'FontWeight', 'normal', 'FontSize', 10);
    
    curLeg = legend(curAx, ...
        'All', 'Onto spines', 'Onto shaft', ...
        'Location', 'EastOutside');
    curLeg.Position([1, 3]) = [0.82, (0.98 - 0.82)];
    curLeg.Box = 'off';
    
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
            'LineWidth', 2, ...
            'FaceAlpha', 1);
	
    curPlot(curAx, curSyns.dist(curSyns.isSpine & ~curSyns.isSoma));
    curPlot(curAx, curSyns.dist(curSyns.isExc));
    
    curAx.TickDir = 'out';
    curAx.Position(3) = 0.8 - curAx.Position(1);
    
    xlim(curAx, curBinEdges([1, end]));
    ylabel(curAx, 'Fraction of synapses');
    ylim(curAx, [0, 1]);
    
    curLeg = legend(curAx, ...
        'Onto spines', ...
        'Excitatory', ...
        'Location', 'EastOutside');
    curLeg.Position([1, 3]) = [0.82, (0.98 - 0.82)];
    curLeg.Box = 'off';
    
    % TC vs CC
    curAx = subplot(4, 1, 3);
    hold(curAx, 'on');
    
    curDiscretize = ...
        @(data) accumarray( ...
            discretize(data, curBinEdges), ...
            1, [numel(curBinEdges) - 1, 1]);
    curBinsAll = curDiscretize(curSyns.dist);
    
    curPlot = ...
        @(ax, data) histogram(ax, ...
            'BinCount', curDiscretize(data), ...
            'BinEdges', curBinEdges, ...
            'Normalization', 'probability', ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2, ...
            'FaceAlpha', 1);
	
    curPlot(curAx, curSyns.dist(curSyns.axonClass == 'Thalamocortical'));
    curPlot(curAx, curSyns.dist(curSyns.axonClass == 'Corticocortical'));
    
    curAx.TickDir = 'out';
    curAx.Position(3) = 0.8 - curAx.Position(1);
    
    xlim(curAx, curBinEdges([1, end]));
    ylabel(curAx, 'Probability');
    
    curLeg = legend(curAx, ...
        'Thalamocortical', ...
        'Corticocortical', ...
        'Location', 'EastOutside');
    curLeg.Position([1, 3]) = [0.82, (0.98 - 0.82)];
    curLeg.Box = 'off';
    
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
        'LineWidth', 2, ...
        'FaceAlpha', 1);
    
    curAx.TickDir = 'out';
    curAx.Position(3) = 0.8 - curAx.Position(1);
    
    ylabel(curAx, 'Median ASI (µm²)');
    xlabel(curAx, 'Distance to soma (µm)');
    xlim(curAx, curBinEdges([1, end]));
    ylim(curAx, [0, 0.4]);
    
    % save figure
    if ~isempty(plotDir)
        curFigFile = sprintf('input-distribution_%s', extWcT.tag{curIdx});
        curFigFile = fullfile(plotDir, curFigFile);
        
        export_fig(strcat(curFigFile, '.png'), curFig);
        export_fig(strcat(curFigFile, '.eps'), curFig);
        clear curFig;
    end
end

%% Orthogonal distance plots
binEdges = linspace(-100, +100, 51);

dimLabels = {'X', 'Y', 'Z'};
synTypes = categories(conn.axonMeta.axonClass);
synTypes = synTypes(1:(end - 1));

for curIdx = size(extWcT, 1)
    curSyns = extWcT.synapses{curIdx};
    if isempty(curSyns); continue; end
    
    curSyns.axonClass = conn.axonMeta.axonClass(curSyns.axonId);
   [~, curSyns.axonClassId] = ismember(curSyns.axonClass, synTypes);
    curSyns(~curSyns.axonClassId, :) = [];
    
    curSyns.dist = extWcT.nodeOrthoDists{curIdx}(curSyns.nodeId, :);
    curSyns.dist = curSyns.dist / 1E3;
    
    curFig = figure();
    curFig.Color = 'white';
    
    for curDimIdx = 1:3
        curData = accumarray( ...
            curSyns.axonClassId, curSyns.dist(:, curDimIdx), ...
            [3, 1], @(dists) {dists}, {zeros(0, 1)});
        curData = cellfun( ...
            @(d) histcounts(d, binEdges), ...
            curData, 'UniformOutput', false);
        curData = transpose(cell2mat(curData));
        
        for curTypeIdx = 1:size(curData, 2)
            curBins = curData(:, curTypeIdx);
            curBinsPerClass = curBins / sum(curBins);
            curBinsPerSyn = curBins / sum(curData(:));
            
            curAx = subplot(3, 3, curDimIdx);
            axis(curAx, 'square');
            hold(curAx, 'on');
            
            histogram(curAx, ...
                'BinEdges', binEdges, ...
                'BinCounts', curBinsPerClass, ...
                'DisplayStyle', 'stairs', ...
                'LineWidth', 2);
            
            curAx = subplot(3, 3, curDimIdx + 3);
            axis(curAx, 'square');
            hold(curAx, 'on');
            
            histogram(curAx, ...
                'BinEdges', binEdges, ...
                'BinCounts', curBinsPerSyn, ...
                'DisplayStyle', 'stairs', ...
                'LineWidth', 2);
        end
        
        curAx.Children = flip(curAx.Children);
        
        curAx = subplot(3, 3, curDimIdx + 2 * 3);
        axis(curAx, 'square');
        hold(curAx, 'on');
        
        curTcEx = curData(:, 2) ./ sum(curData(:, 1:2), 2);
        curTcEx(isnan(curTcEx) | isinf(curTcEx)) = 0;
        curInhEx = curData(:, 3) ./ sum(curData, 2);
        curInhEx(isnan(curInhEx) | isinf(curInhEx)) = 0;
        
        histogram(curAx, ...
            'BinEdges', binEdges, ...
            'BinCounts', curInhEx, ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2);
        histogram(curAx, ...
            'BinEdges', binEdges, ...
            'BinCounts', curTcEx, ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2);
        
        xlabel(curAx, sprintf( ...
            'Synapse location relative to soma (µm along %s)', ...
            dimLabels{curDimIdx}));
    end
    
    curAxes = flip(curFig.Children);
   [curAxes(1:3:end).YLim] = deal([0, max( ...
       arrayfun(@(ax) ax.YLim(end), curAxes(1:3:end)))]);
   [curAxes(2:3:end).YLim] = deal([0, max( ...
       arrayfun(@(ax) ax.YLim(end), curAxes(2:3:end)))]);
   [curAxes(3:3:end).YLim] = deal([0, max( ...
       arrayfun(@(ax) ax.YLim(end), curAxes(3:3:end)))]);
   [curAxes.TickDir] = deal('out');
   
    ylabel(curAxes(1), 'Synapse probability (per axon type)');
    ylabel(curAxes(2), 'Fraction of synaptic input');
    ylabel(curAxes(3), 'Ratio');
   
    curLeg = legend( ...
        curAxes(end - 1), flip(synTypes), ...
        'Location', 'NorthEast');
    curLeg.Box = 'off';
    
    curLeg = legend( ...
        curAxes(end), {'Inh / (Inh + Exc)', 'TC / (TC + CC)'}, ...
        'Location', 'NorthEast');
    curLeg.Box = 'off';

    annotation(curFig, ...
        'textbox', [0, 0.9, 1, 0.1], ...
        'String', { ...
            info.filename; info.git_repos{1}.hash; extWcT.title{curIdx}}, ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center');
end

%% Correlation between ratios and dendrite direction
curMinSyn = 50;
dimLabels = {'X', 'Y', 'Z'};
synTypes = categories(conn.axonMeta.axonClass);

dendT = connectEM.WholeCell.splitWholeCellInputs(wcT, splitNmlT);

dendT.dir = nan(size(dendT.somaPos));
dendT.tcExcRatio = nan(size(dendT.id));
dendT.inhExcRatio = nan(size(dendT.id));

% Let's not analyse dendrites with too few synapses
dendT(cellfun(@height, dendT.synapses) < curMinSyn, :) = [];

for curIdx = 1:size(dendT, 1)
    curSyns = dendT.synapses{curIdx};
    curSomaPos = dendT.somaPos(curIdx, :);
    
    curNodes = dendT.agglo(curIdx).nodes;
    curNodes(isnan(curNodes(:, 4)), :) = [];
    
   [curSegIds, curSegPos] = unique(curNodes(:, 4));
    curSegPos = curNodes(curSegPos, 1:3) - curSomaPos;
    curSegPos = curSegPos .* param.raw.voxelSize;
    
    curSegMass = segMass(curSegIds);
    curSegMass = curSegMass / sum(curSegMass);
    
    curDendDir = curSegPos ./ sqrt(sum(curSegPos .^ 2, 2));
    curDendDir = sum(curSegMass .* curDendDir, 1);
    curDendDir = curDendDir / sqrt(sum(curDendDir .^ 2));
    
    curSynData = accumarray( ...
        double(conn.axonMeta.axonClass(curSyns.axonId)), ...
        1, [numel(synTypes), 1], @sum, 0);
    
    dendT.dir(curIdx, :) = curDendDir;
    dendT.tcExcRatio(curIdx) = curSynData(2) / sum(curSynData(1:2));
    dendT.inhExcRatio(curIdx) = curSynData(3) / sum(curSynData(1:3));
end

curFig = figure();
curFig.Color = 'white';
curFig.Position(3:4) = [690, 510];

for curDimIdx = 1:3
    % TC
    curFit = fit(dendT.dir(:, curDimIdx), dendT.tcExcRatio, 'poly1');
    curAx = subplot(2, 3, curDimIdx);
    
    hold(curAx, 'on');
    scatter(curAx, dendT.dir(:, curDimIdx), dendT.tcExcRatio, 60, '.');
    plot(curAx, [-1, 1], curFit([-1, 1]), 'Color', 'black', 'LineWidth', 2);
    
    % Inh
    curFit = fit(dendT.dir(:, curDimIdx), dendT.inhExcRatio, 'poly1');
    curAx = subplot(2, 3, curDimIdx + 3);
    
    hold(curAx, 'on');
    scatter(curAx, dendT.dir(:, curDimIdx), dendT.inhExcRatio, 60, '.');
    plot(curAx, [-1, 1], curFit([-1, 1]), 'Color', 'black', 'LineWidth', 2);
    
    xlabel(curAx, sprintf('%s-polarity of dendrite', dimLabels{curDimIdx}));
end

curAxes = flip(curFig.Children);
ylabel(curAxes(1), 'TC / (TC + CC)');
ylabel(curAxes(2), 'Inh / (Inh + Exc)');

[curAxes.XLim] = deal([-1, +1]);
[curAxes.PlotBoxAspectRatio] = deal([1, 1, 1]);
[curAxes.DataAspectRatioMode] = deal('auto');

annotation( ...
    curFig, 'textbox', [0, 0.9, 1, 0.1], ...
    'String', {info.filename; info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');

tcExcCorr = arrayfun(@(i) corr(dendT.dir(:, i), dendT.tcExcRatio), 1:3);
tcExcDir = tcExcCorr ./ sqrt(sum(tcExcCorr .^ 2));

inhExcCorr = arrayfun(@(i) corr(dendT.dir(:, i), dendT.inhExcRatio), 1:3);
inhExcDir = inhExcCorr ./ sqrt(sum(inhExcCorr .^ 2));

%% Correlation between inhibition and dendrite length
curMinSyn = 500;
curSynTypes = unique(conn.axonMeta.axonClass);

curWcT = wcT;
curWcT(cellfun(@height, curWcT.synapses) < curMinSyn, :) = [];

curWcT.classConn = cell2mat(cellfun( ...
    @(syn) transpose(accumarray( ...
        conn.axonMeta.axonClassId(syn.axonId), ...
        1, [numel(curSynTypes), 1])), ...
	curWcT.synapses, 'UniformOutput', false));

curWcT.inhFrac = ...
    curWcT.classConn(:, curSynTypes == 'Inhibitory') ...
 ./ sum(curWcT.classConn(:, curSynTypes ~= 'Other'), 2);

curWcT.meanSynDist = cellfun( ...
    @(dists, syn) mean(dists(syn.nodeId)), ...
    curWcT.nodeDists, curWcT.synapses) / 1E3;

curFit = fit(curWcT.meanSynDist, curWcT.inhFrac, 'poly1');

curFig = figure();
curFig.Color = 'white';

curAx = axes(curFig);

hold(curAx, 'on');
scatter(curAx, ...
    curWcT.meanSynDist, curWcT.inhFrac, 60, '.');

curAx.XLim(1) = 0;
curAx.YLim(1) = 0;

plot(curAx, ...
    curAx.XLim, curFit(curAx.XLim), ...
    'Color', 'black', 'LineWidth', 2);

curAx.TickDir = 'out';
xlabel(curAx, 'Mean synapse-to-soma distance (µm)');
ylabel(curAx, 'Inh / (Inh + Exc) ratio');

curTitle = sprintf( ...
	'Whole cells with ≥ %d synapses (n = %d)', ...
    curMinSyn, height(curWcT));
title(curAx, ...
    {info.filename; info.git_repos{1}.hash; curTitle}, ...
    'FontWeight', 'normal', 'FontSize', 10);

%% Quantitative comparison of whole cells
synTypes = categories(conn.axonMeta.axonClass);
wcSynTypes = zeros(size(wcT, 1), numel(synTypes));

for curIdx = 1:size(wcT, 1)
    curSyns = wcT.synapses{curIdx};
    
    curSyns.isSpine = syn.isSpineSyn(curSyns.id);
    curSyns.isSoma = ismember( ...
        wcT.agglo(curIdx).nodes(curSyns.nodeId, 4), ...
        wcT.somaAgglo(curIdx).nodes(:, 4));
    curSyns.type = conn.axonMeta.axonClass(curSyns.axonId);
    
   [~, curSynTypeCount] = ismember(curSyns.type, synTypes);
    curSynTypeCount = accumarray(curSynTypeCount, 1, size(synTypes(:)));
    wcSynTypes(curIdx, :) = curSynTypeCount;
end

%% Plot results
binEdges = linspace(0, 1, 41);

for curWcGroup = wcGroups
    curData = wcSynTypes;
    curData = curData(curWcGroup.wcIds, :);
    
    % Remove cells without synapses
    curData(~any(curData, 2), :) = [];
    curData = curData ./ sum(curData, 2);
    
    curFig = figure();
    curFig.Color = 'white';
    curFig.Position(3:4) = [520, 900];

    for curIdx = 1:numel(synTypes)
        curAx = subplot(numel(synTypes), 1, curIdx);

        histogram( ...
            curAx, curData(:, curIdx), ...
            'BinEdges', binEdges, ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2, ...
            'FaceAlpha', 1);
        
        curAx.XLim = binEdges([1, end]);
        curAx.YLim = [0, size(curData, 1)];
        curAx.TickDir = 'out';
        curAx.Box = 'off';

        title( ...
            curAx, synTypes{curIdx}, ...
            'FontWeight', 'normal', ...
            'FontSize', 10);
    end
    
    curMaxY = [curFig.Children.Children];
    curMaxY = max(arrayfun(@(h) max(h.Values), curMaxY));
   [curFig.Children.YLim] = deal([0, curMaxY]);
   
    xlabel(curAx, 'Fraction of input synapses');
    ylabel(curAx, 'Cells');

    annotation(curFig, ...
        'textbox', [0, 0.9, 1, 0.1], ...
        'String', { ...
            info.filename; info.git_repos{1}.hash; ...
            sprintf('Variability across %s', curWcGroup.title)}, ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center');
end

%% Plot intra-cell variability
binEdges = linspace(0, 1, 21);

synTypes = categories(conn.axonMeta.axonClass);
[~, conn.axonMeta.axonClassId] = ...
    ismember(conn.axonMeta.axonClass, synTypes);

dendT = connectEM.WholeCell.splitWholeCellInputs(wcT, splitNmlT);

extWcT = wcT;
extWcT.dendId(:) = 0;
extWcT = extWcT(:, [1, end, 2:(end - 1)]);
extWcT = cat(1, extWcT, dendT);

curWcIds = unique(extWcT.id);
for curId = reshape(curWcIds, 1, [])
    curDendT = extWcT(extWcT.id == curId, :);
    curDendT = sortrows(curDendT, 'dendId');
    
    curDendT.synData = cell2mat(cellfun(@(s) ...
        transpose(accumarray( ...
            conn.axonMeta.axonClassId(s.axonId), ...
            1, [numel(synTypes), 1])), ...
        curDendT.synapses, 'UniformOutput', false));
    
    curWcT = curDendT(1, :);
    assert(~curWcT.dendId);
    
    curDendT(1, :) = [];
    curDendT(sum(curDendT.synData, 2) < 50, :) = [];
    if size(curDendT, 1) < 2; continue; end
    
    if plotShow
        curFig = figure(); %#ok
    elseif ~isempty(plotDir)
        curFig = figure('visible', 'off');
    else
        % No need to build plot
        continue;
    end
    
    curFig.Color = 'white';
    curFig.Position(3:4) = [720, 410];
    
    curAx = axes(curFig); %#ok
    hold(curAx, 'on');
    
    curData = ...
        curDendT.synData ...
     ./ sum(curDendT.synData, 2);
 
    for curDendIdx = 1:size(curData, 1)
        plot( ...
            curAx, 1:size(curData, 2), curData(curDendIdx, :), ...
            'LineWidth', 2, 'Marker', '.', 'MarkerSize', 24);
    end
    
    plot( ...
        curAx, 1:size(curData, 2), ...
        curWcT.synData ./ sum(curWcT.synData, 2), ...
        'Color', 'black', 'LineStyle', '--', 'LineWidth', 2, ...
        'Marker', '.', 'MarkerSize', 24);
    
    curAx.Box = 'off';
    curAx.TickDir = 'out';
    
    curAx.YLim = [0, 1];
    curAx.XLim = [1 - 0.1, size(curData, 2) + 0.1];
    
    curLeg = cell(size(curDendT, 1) + 1, 1);
    curLeg(1:(end - 1)) = arrayfun( ...
        @(id, nSyn) sprintf('Dendrite %d (n = %d synapses)', id, nSyn), ...
        curDendT.dendId, sum(curDendT.synData, 2), 'UniformOutput', false);
    curLeg{end} = sprintf( ...
        'Full reconstruction (n = %d synapses)', sum(curWcT.synData, 2));
    
    curLeg = legend(curAx, curLeg, 'Location', 'EastOutside');
    curLeg.Box = 'off';
    
    axis(curAx, 'square');
    xticks(curAx, 1:size(curData, 2));
    xticklabels(curAx, synTypes);
    
    xlabel('Synapse type');
    ylabel('Fraction of input synapses');
    
    annotation( ...
        curFig, ...
        'textbox', [0, 0.9, 1, 0.1], ...
    	'String', { ...
            info.filename; info.git_repos{1}.hash; ...
            sprintf('Cell %d', curWcT.id)}, ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    
    if ~isempty(plotDir)
        curFigFile = sprintf('input-variability_%s', curWcT.tag{1});
        curFigFile = fullfile(plotDir, curFigFile);
        
        export_fig(strcat(curFigFile, '.png'), curFig);
        export_fig(strcat(curFigFile, '.eps'), curFig);
        clear curFig;
    end
end

%% Export interesting cells
exportCellIds = [];
exportDir = '/home/amotta/Desktop';

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

for curCellId = reshape(exportCellIds, 1, [])
    curDendT = extWcT(extWcT.id == curCellId & extWcT.dendId ~= 0, :);
    curDendT.agglo = SuperAgglo.connect(curDendT.agglo);
    
    curSkel = ...
        connectEM.WholeCell.inputsOutputsToSkel( ...
            param, curDendT, syn.synapses, segPoints, skel);
        
	curSkelFile = sprintf('cell-%d.nml', curCellId);
    curSkel.write(fullfile(exportDir, curSkelFile));
end
