% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

% Set output directory to write figures to disk instead of displaying them.
plotDir = '';
plotShow = false;

synEmMarginUm = 3;

connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');
wcFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto-and-manual.mat');
somaFile  = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_03_v2.mat');

splitNmlDir = fileparts(fileparts(mfilename('fullpath')));
splitNmlDir = fullfile(splitNmlDir, '+WholeCell', '+Script', 'annotations');

% debugging
debugDir = '';
debugCellIds = [];

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

segMass = Seg.Global.getSegToSizeMap(param);
segPoint = Seg.Global.getSegToPointMap(param);

[conn, syn, axonClasses] = connectEM.Connectome.load(param, connFile);

synPos = connectEM.Synapse.calculatePositions(param, syn);
synTypes = unique(conn.axonMeta.axonClass);

wcData = load(wcFile);
somaData = load(somaFile);

%% NML files for splitting whole cells into individual dendrites
splitNmlT = connectEM.WholeCell.loadSplitNmls(splitNmlDir);
splitNmlT.cellId = wcData.idxWholeCells(splitNmlT.aggloId);
assert(all(splitNmlT.cellId));

%% Complete whole cell somata
clear cur*;

wcT = table;
wcT.id = wcData.idxWholeCells(wcData.indWholeCells);
wcT.agglo = wcData.dendrites(wcData.indWholeCells);
SuperAgglo.check(wcT.agglo);

% Find corresponding somata
[~, curSomaIds] = ismember(wcT.id, somaData.idxSomata);
wcT.somaAgglo = somaData.dendrites(curSomaIds);

curCalcSomaPos = @(n) ...
    sum(segMass(n(:, 4)) .* n(:, 1:3), 1) ....
    ./ sum(segMass(n(:, 4)));
wcT.somaPos = cell2mat(arrayfun( ...
    @(a) curCalcSomaPos(a.nodes), ...
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

%% Calculate dendritic path length to soma
% Nodes within the soma are assigned as distance of zero.
clear cur*;

wcT.nodeDists = cell(size(wcT.id));

for curIdx = 1:height(wcT)
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
    
    % NOTE(amotta): Travelling within soma doesn't add path length.
    curDists(all(curSomaNodeMask(curAgglo.edges), 2)) = 0;
    
    curDists = graph( ...
        curAgglo.edges(:, 1), ...
        curAgglo.edges(:, 2), ...
        curDists, curNodeCount);
    curDists = distances(curDists, curSomaNodeId);
    
    % Sanity checks
    assert(~any(curDists(curSomaNodeMask)));
    assert(all(curDists(~curSomaNodeMask)));
    
    curDists = reshape(curDists, [], 1);
    wcT.nodeDists{curIdx} = curDists;
end

%% Collect input synapses
clear cur*;

wcT.synapses = cell(size(wcT.id));
wcT.classConn = nan(height(wcT), numel(synTypes));

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
    
    curSynT.nodeId = cellfun( ...
        @(ids) intersect(ids, curAgglo.nodes(:, 4)), ...
        syn.synapses.postsynId(curSynT.id), ...
        'UniformOutput', false);
    
    % Assign synapse to largest segment
   [~, curNodeId] = cellfun( ...
        @(ids) max(segMass(ids)), curSynT.nodeId);
    curNodeIds = arrayfun( ...
        @(ids, idx) ids{1}(idx), curSynT.nodeId, curNodeId);
    
    % Translate segment to node IDs
   [~, curSynT.nodeId] = ismember( ...
       double(curNodeIds), curAgglo.nodes(:, 4));
    
    curSynT.axonId = repelem( ...
        curConnRows.edges(:, 1), ...
        cellfun(@numel, curConnRows.synIdx));
    curClassConn = accumarray( ...
        double(conn.axonMeta.axonClass(curSynT.axonId)), ...
        1, [numel(synTypes), 1]);
    
    wcT.synapses{curIdx} = curSynT;
    wcT.classConn(curIdx, :) = curClassConn;
end

%% Debugging
clear cur*;

if ~isempty(debugDir)
    skel = skeleton();
    skel = Skeleton.setParams4Pipeline(skel, param);
    skel = skel.setDescription(sprintf( ...
        '%s (%s)', info.filename, info.git_repos{1}.hash));

    for curIdx = 1:numel(debugCellIds)
        curId = debugCellIds(curIdx);
        curAgglo = wcT.agglo(curId);
        
        curSynDists = wcT.nodeDists{curId} / 1E3;
        curSynDists = curSynDists(wcT.synapses{curId}.nodeId);
        curSynPos = synPos(wcT.synapses{curId}.id, :);

        % generate skeleton
        curSkel = skel.addTree( ...
            sprintf('Whole cell %d', curId), ...
            curAgglo.nodes(:, 1:3), curAgglo.edges);
        curSkel = curSkel.addNodesAsTrees(curSynPos);
        
        % add distance annotations
        curSkel.names(2:end) = arrayfun( ...
            @(distUm) sprintf('Synapse. %.1f µm to soma.', distUm), ...
            curSynDists, 'UniformOutput', false);

        curSkelName = sprintf( ...
            '%0*d_whole-cell-%d.nml', ...
            ceil(log10(1 + size(wcT, 1))), curIdx, curId);
        curSkel.write(fullfile(debugDir, curSkelName));
    end
end

%% Remove interneurons and exc. synapses onto exc. somata
% As discussed with MH on 08.08.2018

% Remove interneurons
interNeuronCellIds = conn.denMeta.cellId(conn.denMeta.isInterneuron);
interNeuronCellIds = setdiff(interNeuronCellIds, 0);

wcT(ismember(wcT.id, interNeuronCellIds), :) = [];
splitNmlT(ismember(splitNmlT.cellId, interNeuronCellIds), :) = [];

% Remove excitatory synapses onto somata
for curIdx = 1:height(wcT)
    curSynapses = wcT.synapses{curIdx};
    curNodeDists = wcT.nodeDists{curIdx};
    
    curExcMask = ismember( ...
        conn.axonMeta.axonClass(curSynapses.axonId), ...
        {'Corticocortical', 'Thalamocortical'});
    curSomaMask = ~curNodeDists(curSynapses.nodeId);
    curDropMask = curExcMask & curSomaMask;
    curSynapses(curDropMask, :) = [];
    
    % Update class connectome
    curClassConn = double(conn.axonMeta.axonClass(curSynapses.axonId));
    curClassConn = accumarray(curClassConn, 1, [numel(synTypes), 1]);
    
    wcT.synapses{curIdx} = curSynapses;
    wcT.classConn(curIdx, :) = curClassConn;
end

%% Render isosurface
% Since Mr. Amira is on vacation.
clear cur*;

% Configuration
curIsoDir = '/tmpscratch/amotta/l4/2018-05-10-whole-cell-isosurfaces-spine-evolution/full/mat/';

curCellId = 21;
curCamPos = [44683, 16045, -15604];
curCamTarget = [2852.5, 4320.0, 1965.4];
curCamUpVector = [-0.0178, -0.0785, -0.0153];
curCamViewAngle = 10.4142;

curIso = load(fullfile(curIsoDir, sprintf('iso-%d.mat', curCellId)));
curIso = curIso.isoSurf;

curSyn = wcT.synapses{wcT.id == curCellId};
curSyn.type = conn.axonMeta.axonClass(curSyn.axonId);

curSynPos = connectEM.Synapse.calculatePositions(param, syn, 'pre');
curSynPos = curSynPos(curSyn.id, :);

curFig = figure;
curFig.Color = 'none';

curAx = axes(curFig);
curAx.Visible = 'off';

curColors = curAx.ColorOrder;
curColors = cat(1, curColors(1:3, :), [0, 0, 0]);
curRads = [1, 1, 1, 0.5] .* 1E3;

hold(curAx, 'on');
daspect(curAx, 1 ./ param.raw.voxelSize);

curP = patch(curAx, curIso);
curP.EdgeColor = 'none';
curP.FaceColor = repelem(0.85, 3);
material(curP, 'dull');

[curX, curY, curZ] = sphere();
curX = curX / param.raw.voxelSize(1);
curY = curY / param.raw.voxelSize(2);
curZ = curZ / param.raw.voxelSize(3);

for curId = 1:height(curSyn)
    curTypeId = double(curSyn.type(curId));
    curColor = curColors(curTypeId, :);
    curPos = curSynPos(curId, :);
    curRad = curRads(curTypeId);
    
    curSurf = surf(curAx, ...
        curRad * curX + curPos(1), ...
        curRad * curY + curPos(2), ...
        curRad * curZ + curPos(3));
    curSurf.FaceColor = curColor;
    curSurf.EdgeColor = 'none';
    material(curSurf, 'dull');
end

curAx.CameraPosition = curCamPos;
curAx.CameraTarget = curCamTarget;
curAx.CameraUpVector = curCamUpVector;
curAx.CameraViewAngle = curCamViewAngle;

camlight(curAx);

curScaleBar = 10E3;
curScaleBar = curScaleBar ./ param.raw.voxelSize;

curOrig = [ ...
    param.bbox(1, 1), ...
    param.bbox(2, 2) - 1E3, ...
    param.bbox(3, 1) + 3E3];
curDirs = [+1, -1, +1];

plot3(curAx, ...
    curOrig(1) + curDirs(1) .* [0, curScaleBar(1)], ...
    repelem(curOrig(2), 2), ...
    repelem(curOrig(3), 2));
plot3(curAx, ...
    repelem(curOrig(1), 2), ...
    curOrig(2) + curDirs(2) .* [0, curScaleBar(2)], ...
    repelem(curOrig(3), 2), 'Color', 'black');
plot3(curAx, ...
    repelem(curOrig(1), 2), ...
    repelem(curOrig(2), 2), ...
    curOrig(3) + curDirs(3) .* [0, curScaleBar(3)]);

curPlots = findobj(curAx, 'Type', 'Line');
set(curPlots, 'Color', 'black', 'LineWidth', 2);

annotation( ...
    curFig, 'textbox', [0, 0.9, 1, 0.1], ...
    'String', {info.filename; info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');

%% Plot soma position versus input ratios
clear cur*;
curMinSyn = 100;
curBinSizeUm = 2;

curLimits = (param.bbox(1, 2) - param.bbox(1, 1)) / 2;
curLimits = curLimits * param.raw.voxelSize(1) / 1E3;
curLimY = floor(curLimits) .* [-1, +1];

curBinEdges = curLimits - synEmMarginUm;
curBinEdges = curBinSizeUm * floor(curBinEdges / curBinSizeUm);
curBinEdges = (-curBinEdges):curBinSizeUm:curBinEdges;

curSynT = connectEM.Connectome.buildSynapseTable(conn, syn);
curSynT.synType = conn.axonMeta.axonClass(curSynT.preAggloId);

curSynPos = synPos;
curSynPos = curSynPos(:, 1) - mean(param.bbox(1, :));
curSynPos = curSynPos * param.raw.voxelSize(1) / 1E3;

curSynT.pos = curSynPos(curSynT.id);
curSynT.binId = discretize(curSynT.pos, curBinEdges);
curSynT(isnan(curSynT.binId), :) = [];

curBinCounts = accumarray( ...
    cat(2, curSynT.binId, double(curSynT.synType)), ...
    1, [numel(curBinEdges) - 1, numel(synTypes)]);
curBinRatios = curBinCounts(:, 2) ./ sum(curBinCounts(:, 1:2), 2);

curFitRatio = (curBinEdges(1:(end - 1)) + curBinEdges(2:end)) / 2;
curFitRatio = fit(curFitRatio(:), curBinRatios(:), 'poly1');

curWcT = wcT;
curWcT(sum(curWcT.classConn, 2) < curMinSyn, :) = [];

curWcT.tcRatio = ...
    curWcT.classConn(:, 2) ...
 ./ sum(curWcT.classConn(:, 1:2), 2);

curWcT.somaPosRel = curWcT.somaPos - transpose(mean(param.bbox, 2));
curWcT.somaPosRel = curWcT.somaPosRel .* param.raw.voxelSize / 1E3;

% Correct for YZ changes
curFit = curWcT.tcRatio - mean(curWcT.tcRatio);
curFit = fit(curWcT.somaPosRel(:, 2:3), curFit, 'poly11');
curWcT.corrTcRatio = curWcT.tcRatio - curFit(curWcT.somaPosRel(:, 2:3));

curFig = figure();
curFig.Color = 'white';
curFig.Position(3:4) = [290, 200];

% Correct for YZ
curAx = axes(curFig);
hold(curAx, 'on');

curAx.ColorOrder = [ ...
    0.4660, 0.6740, 0.1880;
    0.0000, 0.4470, 0.7410];

scatter(curAx, curWcT.corrTcRatio, curWcT.somaPosRel(:, 1), 60, '.');
scatter(curAx, curWcT.tcRatio, curWcT.somaPosRel(:, 1), 60, '.');

curFit = fit( ...
    curWcT.somaPosRel(:, 1), ...
    curWcT.corrTcRatio, 'poly1');
plot(curAx, ...
    curFitRatio(curLimY), curLimY, ...
    'Color', 'red', 'LineWidth', 2, 'LineStyle', '--');
plot(curAx, ...
    curFit(curLimY), curLimY, ...
    'Color', 'black', 'LineWidth', 2);

xlabel(curAx, {'TC / (TC + CC)'; 'per neuron'});

curAx.TickDir = 'out';
curAx.YDir = 'reverse';
curAx.YLim = curLimY;

curLeg = {'YZ corrected ratios'; 'Raw ratios'};
curLeg = legend(curAx, curLeg, 'Location', 'EastOutside');
curLeg.Box = 'off';

title(curAx, ...
    {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);

%% Correlation between ratios and dendrite direction
curMinSyn = 50;

dendT = connectEM.WholeCell.splitWholeCellInputs(wcT, splitNmlT);
dendT(cellfun(@height, dendT.synapses) < curMinSyn, :) = [];

[~, dendT.cellRow] = ismember(dendT.id, curWcT.id);
dendT(~dendT.cellRow, :) = [];

dendT.dir = nan(size(dendT.somaPos));
dendT.classConn = nan(height(dendT), numel(synTypes));
dendT.tcExcRatio = nan(size(dendT.id));
dendT.wcRelTcExcRatio = nan(size(dendT.id));

for curIdx = 1:size(dendT, 1)
    curSyns = dendT.synapses{curIdx};
    curCell = curWcT(dendT.cellRow(curIdx), :);
    
    % Calculate dendrite orientation
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
    dendT.classConn(curIdx, :) = curSynData;
    
    dendT.tcExcRatio(curIdx) = ...
        curSynData(2) / sum(curSynData(1:2));
    dendT.wcRelTcExcRatio(curIdx) = ...
        dendT.tcExcRatio(curIdx) / curCell.tcRatio;
end

%% Plotting
curVars = { ...
    'wcRelTcExcRatio', {'TC / (TC + CC)'; 'normalized to cell'}};
curVarLabels = curVars(:, 2);
curVars = curVars(:, 1);

curBinEdges = [-1, -0.5, 0.5, 1];

curFig = figure();
curFig.Color = 'white';
curFig.Position(3:4) = [400, 220];

for curVarIdx = 1:numel(curVars)
    curVar = curVars{curVarIdx};
    curVarLabel = curVarLabels{curVarIdx};
    
    curFracs = dendT.(curVar);
    curDirs = dendT.dir(:, 1);
    curBinIds = discretize(curDirs, curBinEdges);
    
    % Statistical tests
    fitlm(curDirs, curFracs) %#ok
   [~, curPVal] = ttest2( ...
       curFracs(curBinIds == 1), ...
       curFracs(curBinIds == 3)) %#ok
    
    curFit = fit(curDirs, curFracs, 'poly1');
    
    curAx = subplot(numel(curVars), 2, (curVarIdx - 1) * 2 + 1);
    hold(curAx, 'on');
    
    scatter(curAx, curFracs, curDirs, '.');
    plot(curAx, ...
        curFit([-1, +1]), [-1, +1], ...
        'Color', 'black', 'LineWidth', 2);
    
    xlabel(curAx, curVarLabel);
    ylabel(curAx, 'Alignment to cortical axis');
    curAx.YTickLabel{1} = sprintf('%s', curAx.YTickLabel{1});
    curAx.YTickLabel{end} = sprintf('%s', curAx.YTickLabel{end});
    
    curAx = subplot(numel(curVars), 2, (curVarIdx - 1) * 2 + 2);
    boxplot( ...
        curAx, curFracs, curBinIds, ...
        'Orientation', 'horizontal', 'Symbol', '.', ...
        'Positions', (curBinEdges(1:(end - 1)) + curBinEdges(2:end)) / 2, ...
        'Widths', (curBinEdges(2:end) - curBinEdges(1:(end - 1))) / 2);
    xlabel(curAx, curVarLabel);
    ylabel(curAx, 'Alignment to cortical axis');
    
    curAx.Box = 'off';
    curAx.YLim = [-1, +1];
    yticks(curAx, curBinEdges);
    yticklabels(curAx, arrayfun(@num2str, ...
        curBinEdges, 'UniformOutput', false));
end

curAx = flip(cat(1, curFig.Children));
set(curAx, ...
    'XLim', [0, 2], ...
    'YDir', 'reverse', ...
    'TickDir', 'out', ...
    'PlotBoxAspectRatio', [1, 1, 1], ...
    'DataAspectRatioMode', 'auto');

annotation( ...
    curFig, 'textbox', [0, 0.9, 1, 0.1], ...
    'String', {info.filename; info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');

% MATLAB's `boxplot` function is a nightmare
curAx(2).Position(2:4) = curAx(1).Position(2:4);

%% Statistical significance test
fprintf([ ...
    'Significance test for dendrite / neuron ', ...
    'tc / (tc + cc) ratio\n']);
[~, curTopBotPVal] = ttest2( ...
    dendT.wcRelTcExcRatio(dendT.dir(:, 1) < -0.5), ...
    dendT.wcRelTcExcRatio(dendT.dir(:, 1) > +0.5));
fprintf([ ...
    '→ Pia- / WM-oriented dendrites different ', ...
    'with p-value of %g\n'], curTopBotPVal)
curFit = fitlm(dendT.dir(:, 1), dendT.wcRelTcExcRatio);
disp(curFit);
fprintf('\n');
fprintf('\n');

tcCcCellNormalizedWmPiaRatio = ...
    curFit.predict(+1) / curFit.predict(-1) %#ok

curFit = fitlm(dendT.dir(:, 1), dendT.tcExcRatio);
tcCcWmPiaRatio = ...
    curFit.predict(+1) / curFit.predict(-1) %#ok

%% Try to find border cells
clear cur*;
somaPos = cell2mat(arrayfun( ...
    @(s) mean(s.nodes(:, 1:3), 1), ...
    wcT.somaAgglo, 'UniformOutput', false));

voxelSize = param.raw.voxelSize;
curSomaDistXY = voxelSize(3) * min(pdist2( ...
    somaPos(:, 3), transpose(param.bbox(3, :))), [], 2);
curSomaDistXZ = voxelSize(2) * min(pdist2( ...
    somaPos(:, 2), transpose(param.bbox(2, :))), [], 2);
curSomaDistYZ = voxelSize(1) * min(pdist2( ...
    somaPos(:, 1), transpose(param.bbox(1, :))), [], 2);

curYzWcIds = find(curSomaDistYZ >= 10E3 & ...
    (curSomaDistXY < 10E3 | curSomaDistXZ < 10E3));

wcGroups = struct;
wcGroups(1).wcRows = curYzWcIds;
wcGroups(1).title = sprintf( ...
    'whole cells in YZ view (n = %d)', ...
    numel(wcGroups(1).wcRows));
wcGroups(1).tag = 'yz-view';

wcGroups(2).wcRows = find( ...
    cellfun(@height, wcT.synapses) > 100);
wcGroups(2).title = sprintf( ...
    'whole cells with > 100 synapses (n = %d)', ...
    numel(wcGroups(2).wcRows));
wcGroups(2).tag = '100-syn';

wcGroups(3).wcRows = find( ...
    cellfun(@height, wcT.synapses) > 0);
wcGroups(3).title = sprintf( ...
    'whole cells with > 0 synapses (n = %d)', ...
    numel(wcGroups(3).wcRows));
wcGroups(3).tag = '1-syn';

%% Generate queen neurons
extWcT = wcT;
extWcT.nodeParentIds = arrayfun( ...
    @(id, nodeCount) reshape(repelem(id, nodeCount), [], 1), ...
    extWcT.id, cellfun(@numel, extWcT.nodeDists), ...
    'UniformOutput', false);

for curIdx = 1:numel(wcGroups)
    curWcGroup = wcGroups(curIdx);
    curWcT = extWcT(curWcGroup.wcRows, :);
    
    curNodeIdOff = cumsum(cellfun(@numel, curWcT.nodeDists));
    curNodeIdOff = [0; curNodeIdOff(1:(end - 1))];
    
    curLumpedSynapses = cat(1, curWcT.synapses{:});
    curLumpedSynapses.nodeId = curLumpedSynapses.nodeId ...
        + repelem(curNodeIdOff, cellfun(@height, curWcT.synapses));

    curQueenWcT = table;
    curQueenWcT.id = 0;

    curQueenWcT.agglo = struct;
    curQueenWcT.agglo.nodes = cat(1, curWcT.agglo.nodes);
    curQueenWcT.agglo.edges = zeros(0, 2);

    curQueenWcT.somaAgglo = struct;
    curQueenWcT.somaAgglo.nodes = cat(1, curWcT.somaAgglo.nodes);
    curQueenWcT.somaAgglo.edges = zeros(0, 2);
    
    curQueenWcT.somaPos = nan(1, 3);
    
    curQueenWcT.title = curWcGroup.title;
    curQueenWcT.tag = curWcGroup.tag;

    curQueenWcT.nodeDists = {cell2mat(curWcT.nodeDists)};
    curQueenWcT.nodeParentIds = {cat(1, curWcT.nodeParentIds{:})};
    curQueenWcT.synapses = {curLumpedSynapses};
    curQueenWcT.classConn = sum(curWcT.classConn, 1);

    extWcT = cat(1, extWcT, curQueenWcT);
end

%% Numbers
curDists = [40, 70];
curBinSize = 10;

curIdx = height(extWcT);
curSyns = extWcT.synapses{curIdx};
curSyns.axonClass = conn.axonMeta.axonClass(curSyns.axonId);
curSyns.dist = extWcT.nodeDists{curIdx}(curSyns.nodeId);
curSyns.dist = curSyns.dist / 1E3;

for curDist = curDists
    curBinEdges = curDist + curBinSize .* [-1, +1] ./ 2;
    curSynMask = ...
        curSyns.dist >= curBinEdges(1) ...
      & curSyns.dist <= curBinEdges(2);

    curSynTypeData = accumarray( ...
        double(curSyns.axonClass(curSynMask)), 1, [4, 1]);
    curInhExcRatio = curSynTypeData(3) / sum(curSynTypeData(1:3));

    fprintf('Distance from soma: %g µm\n', curDist);
    fprintf('* Inh / (Inh + Exc): %g\n', curInhExcRatio);
end

%% Plotting
globalRatios = connectEM.Connectome.buildSynapseTable(conn, syn);
globalRatios.axonClass = conn.axonMeta.axonClass(globalRatios.preAggloId);
globalRatios = accumarray(double(globalRatios.axonClass), 1, [4, 1]);

globalRatios = [ ...
    globalRatios(3) / sum(globalRatios(1:3)), ...
    globalRatios(2) / sum(globalRatios(1:2))];

for curIdx = height(extWcT)
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
    curFig.Position(3:4) = [220, 220];
    
    curSyns.parentId = extWcT.nodeParentIds{curIdx}(curSyns.nodeId);
    curSyns.dist = extWcT.nodeDists{curIdx}(curSyns.nodeId);
    curSyns.dist = curSyns.dist / 1E3;
    
    % Move soma synapses to separate bin
    curSyns.isSoma = ismember( ...
        extWcT.agglo(curIdx).nodes(curSyns.nodeId, 4), ...
        extWcT.somaAgglo(curIdx).nodes(:, 4));
    curSyns.dist(curSyns.isSoma) = -eps;
    
    curMaxDist = 10 * ceil(max(curSyns.dist) / 10);
    curBinEdges = -10:10:curMaxDist;
    
    curPlotRange = 10 * ceil(prctile(curSyns.dist, 99) / 10);
    curPlotRange = [curBinEdges(1), curPlotRange]; %#ok
    
    curSyns.axonClass = conn.axonMeta.axonClass(curSyns.axonId);
    
    curBinIds = discretize(curSyns.dist, curBinEdges);
    curSynTypeData = accumarray( ...
        cat(2, curBinIds, double(curSyns.axonClass)), ...
        1, [numel(curBinEdges) - 1, 4]);
    curNeuronCounts = accumarray( ...
        curBinIds, curSyns.parentId, ...
        [numel(curBinEdges) - 1, 1], ...
        @(ids) numel(unique(ids)));
    
    curPlot = @(ax, varargin) ...
        histogram( ...
            ax, varargin{:}, ...
            'BinEdges', curBinEdges, ...
            'DisplayStyle', 'stairs', ...
            'LineWidth', 2, ...
            'FaceAlpha', 1);
    
    curAx = subplot(3, 1, 1);
    hold(curAx, 'on');
        
	curPlot(curAx, curSyns.dist(curSyns.axonClass == 'Corticocortical'));
	curPlot(curAx, curSyns.dist(curSyns.axonClass == 'Thalamocortical'));
	curPlot(curAx, curSyns.dist(curSyns.axonClass == 'Inhibitory'));
    
    curAx.TickDir = 'out';
    curAx.Position(3) = 0.8 - curAx.Position(1);
    
    xlim(curAx, curPlotRange);
    ylabel(curAx, 'Synapses');
    
    curLeg = legend(curAx, ...
        'Corticocortical', ...
        'Thalamocortical', ...
        'Inhibitory', ...
        'Location', 'EastOutside');
    curLeg.Position([1, 3]) = [0.82, (0.98 - 0.82)];
    curLeg.Box = 'off';
    
    title(curAx, { ...
        sprintf('Inputs onto %s', extWcT.title{curIdx}); ...
        sprintf('%s (%s)', info.filename, info.git_repos{1}.hash)}, ...
        'FontWeight', 'normal', 'FontSize', 10);
        
    curAx = subplot(3, 1, 2);
    hold(curAx, 'on');
    
    curPlot(curAx, 'BinCounts', curNeuronCounts, 'EdgeColor', 'black');
    
    curAx.YDir = 'reverse';
    curAx.Position(3) = 0.8 - curAx.Position(1);
    curAx.TickDir = 'out';
    
    xlim(curAx, curPlotRange);
    ylabel(curAx, 'Neurons');
    
    % Synapse ratios
    curAx = subplot(3, 1, 3);
    hold(curAx, 'on');
	
    curTemp = curSynTypeData(:, 3) ./ sum(curSynTypeData(:, 1:3), 2);
    curTemp(isnan(curTemp)) = 0;
    curPlot(curAx, 'BinCount', curTemp);
    
    curTemp = curSynTypeData(:, 2) ./ sum(curSynTypeData(:, 1:2), 2);
    curTemp(isnan(curTemp)) = 0;
    curPlot(curAx, 'BinCount', curTemp);
    
    curAx.TickDir = 'out';
    curAx.Position(3) = 0.8 - curAx.Position(1);
    
    ylim(curAx, [0, 1]);
    xlim(curAx, curPlotRange);
    ylabel(curAx, 'Ratio');
    
    plot(curAx, ...
        curAx.XLim, repmat(globalRatios(1), 1, 2), ...
        'Color', curAx.ColorOrder(1, :), 'LineStyle', '--');
    plot(curAx, ...
        curAx.XLim, repmat(globalRatios(2), 1, 2), ...
        'Color', curAx.ColorOrder(2, :), 'LineStyle', '--');
    
    curLeg = legend(curAx, ...
        'Inh / (Inh + Exc)', ...
        'TC / (TC + CC)', ...
        'Location', 'EastOutside');
    curLeg.Position([1, 3]) = [0.82, (0.98 - 0.82)];
    curLeg.Box = 'off';
    
    curAxes = findobj(curFig, 'Type', 'Axes');
    
    curTickLabels = repelem({''}, numel(curBinEdges));
    curTickLabels{curBinEdges == 0} = '0';
    curTickLabels{end} = num2str(curBinEdges(end));
    
    set(curAxes, ...
        'Color', 'none', ...
        'XLim', curBinEdges([1, end]), ...
        'XTick', curBinEdges, ...
        'XTickLabels', curTickLabels);
    
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
            curAx, char(synTypes(curIdx)), ...
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
            double(conn.axonMeta.axonClass(s.axonId)), ...
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
    xticklabels(curAx, arrayfun( ...
        @char, synTypes, 'UniformOutput', false));
    
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
