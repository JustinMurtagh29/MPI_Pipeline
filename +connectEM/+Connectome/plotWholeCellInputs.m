% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

% Set output directory to write figures to disk instead of displaying them.
plotDir = '/home/amotta/Desktop/whole-cell-inputs';
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

[conn, syn] = connectEM.Connectome.load(param, connFile);

wcData = load(wcFile);
somaData = load(somaFile);

%% NML files for whole cell splitting
splitNmlT = table;
splitNmlT.path = dir(fullfile(splitNmlDir, '*.nml'));
splitNmlT.path = reshape({splitNmlT.path.name}, [], 1);
splitNmlT.path = fullfile(splitNmlDir, splitNmlT.path);

splitNmlT.nml = cellfun(@slurpNml, splitNmlT.path);

splitNmlT.cellId = arrayfun( ...
    @(nml) regexp( ...
        nml.parameters.experiment.description, ...
        'Agglomerate (\d+)$', 'tokens', 'once'), ...
    splitNmlT.nml);
splitNmlT.cellId = str2double(splitNmlT.cellId);
splitNmlT.cellId = wcData.idxWholeCells(splitNmlT.cellId);
assert(all(splitNmlT.cellId));

splitNmlT.dendNodes = cell(size(splitNmlT.path));
for curIdx = 1:numel(splitNmlT.dendNodes)
    curNml = splitNmlT.nml(curIdx);
    curTrees = NML.buildTreeTable(curNml);
    curNodes = NML.buildNodeTable(curNml);
    curComments = NML.buildCommentTable(curNml);
    
    curNodes.nodeId(:) = {''};
   [curMask, curIds] = ismember(curNodes.id, curComments.node);
    curNodes.nodeId(curMask) = curComments.comment(curIds(curMask));
    
    curNodes.nodeId = cellfun( ...
        @(nodeId) regexp( ...
            nodeId, 'Node (\d+)$', 'tokens', 'once'), ...
        curNodes.nodeId, 'UniformOutput', false);
    curNodes.nodeId(cellfun(@isempty, curNodes.nodeId)) = {{'0'}};
    curNodes.nodeId = str2double(vertcat(curNodes.nodeId{:}));
    curNodes(~curNodes.nodeId, :) = [];
    
    % Remove soma tree(s)
    curSomaTreeIds = curTrees.id(contains( ...
        curTrees.name, 'soma', 'IgnoreCase', true));
    curNodes(ismember(curNodes.treeId, curSomaTreeIds), :) = [];
    
   [~, ~, curNodes.dendId] = unique(curNodes.treeId);
    splitNmlT.dendNodes{curIdx} = accumarray( ...
        curNodes.dendId, curNodes.nodeId, ...
        [max([1; curNodes.dendId(:)]), 1], ...
        @(ids) {ids}, {zeros(0, 1)});
end

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

    curQueenWcT.nodeDists = {cell2mat(curWcT.nodeDists)};
    curQueenWcT.synapses = {lumpedSynapses};
    
    curQueenWcT.title = curWcGroup.title;
    curQueenWcT.tag = curWcGroup.tag;

    extWcT = cat(1, extWcT, curQueenWcT);
end

%% Plotting
for curIdx = 1:size(extWcT, 1)
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
binEdges = linspace(0, 1, 21);

for curWcGroup = wcGroups
    curData = wcSynTypes;
    curData = curData(curWcGroup.wcIds, :);
    
    % Remove cells without synapses
    curData(~any(curData, 2), :) = [];
    curData = curData ./ sum(curData, 2);

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
