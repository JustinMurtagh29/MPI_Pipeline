% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
isoDir = '/tmpscratch/amotta/l4/2018-05-10-whole-cell-isosurfaces-spine-evolution/full/mat/';

% Set output directory to write figures to disk instead of displaying them.
plotDir = '';
plotShow = true;

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

synTypeFracs = connectEM.Connectome.buildSynapseTable(conn, syn);
synTypeFracs = conn.axonMeta.axonClass(synTypeFracs.preAggloId);
[~, synTypeFracs] = ismember(synTypeFracs, synTypes);
synTypeFracs = accumarray(synTypeFracs, 1);
synTypeFracs = synTypeFracs / sum(synTypeFracs);

wcData = load(wcFile);
somaData = load(somaFile);

colors = get(groot, 'defaultAxesColorOrder');

%% Find axon subpopulations with target specificities
clear cur*;
specConn = conn;
specAxonClasses = axonClasses(1:2);

[specConn, specAxonClasses] = ...
    connectEM.Connectome.prepareForSpecificityAnalysis( ...
        specConn, specAxonClasses, 'minSynPre', 10);
specAxonClasses = ...
    connectEM.Connectome.buildAxonSpecificityClasses( ...
        specConn, specAxonClasses);

% Reset classes in axon meta data. The values are filled in below.
specSynTypes = {};
specConn.axonMeta.axonClass = ...
    repelem({[]}, numel(specConn.axons), 1);

% Modified from
% +connectEM/+Consistency/loadConnectome.m
% commit a23788a05313f61dfef152443cde5d7ff59d7622
for curIdx = 1:numel(specAxonClasses)
    curAxonClass = specAxonClasses(curIdx);
    curSpecs = curAxonClass.specs;
    
    curPrefix = strsplit(curAxonClass.title);
    curPrefix = curPrefix{1};
    curPrefix(1) = upper(curPrefix(1));
    
    curNames = fieldnames(curSpecs);
    curAxonIds = cellfun( ...
        @(n) curSpecs.(n).axonIds, ...
        curNames, 'UniformOutput', false);

    % Find axons with no and multiple specificities
   [curSpecAxonIds, ~, curDupAxonIds] = unique(cell2mat(curAxonIds));
    curDupAxonIds = curSpecAxonIds(accumarray(curDupAxonIds, 1) > 1);
    curNonSpecAxonIds = setdiff(curAxonClass.axonIds, curSpecAxonIds);

    curAxonIds = cellfun( ...
        @(ids) setdiff(ids, curDupAxonIds), ...
        curAxonIds, 'UniformOutput', false);
    
    curNames{end + 1} = 'MultiSpec'; %#ok
    curAxonIds{end + 1} = curDupAxonIds(:); %#ok
    
    curNames{end + 1} = 'NoSpec'; %#ok
    curAxonIds{end + 1} = curNonSpecAxonIds(:); %#ok
    
    for curNameIdx = 1:numel(curNames)
        curName = curNames{curNameIdx};
        curSpecs.(curName).axonIds = curAxonIds{curNameIdx};
    end
    
    % Update specificity classes
    specAxonClasses(curIdx).specs = curSpecs;
    
    % Update axon meta data
    curNames = strcat(curPrefix, curNames);
    specSynTypes = cat(1, specSynTypes(:), curNames(:));
    specConn.axonMeta.axonClass(cell2mat(curAxonIds)) = ...
        repelem(curNames, cellfun(@numel, curAxonIds), 1);
end

specSynTypes{end + 1} = 'Other';
specConn.axonMeta.axonClass(cellfun( ...
    @isempty, specConn.axonMeta.axonClass)) = {'Other'};

specSynTypes = categorical(specSynTypes, specSynTypes, 'Ordinal', true);
[~, curIds] = ismember(specConn.axonMeta.axonClass, specSynTypes);
specConn.axonMeta.axonClass = specSynTypes(curIds);

curSynT = connectEM.Connectome.buildSynapseTable(specConn, syn);
curAxonIds = conn.axonMeta.axonClass(curSynT.preAggloId);
[~, curAxonIds] = ismember(curAxonIds, synTypes);
curSpecIds = specConn.axonMeta.axonClass(curSynT.preAggloId);
[~, curSpecIds] = ismember(curSpecIds, specSynTypes);

specAxonTypeFracs = accumarray([curAxonIds, curSpecIds], 1);
specAxonTypeFracs = specAxonTypeFracs / sum(specAxonTypeFracs(:));

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
wcT.specClassConn = nan(height(wcT), numel(specSynTypes));

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
    curSpecClassConn = accumarray( ...
        double(specConn.axonMeta.axonClass(curSynT.axonId)), ...
        1, [numel(specSynTypes), 1]);
    
    wcT.synapses{curIdx} = curSynT;
    wcT.classConn(curIdx, :) = curClassConn;
    wcT.specClassConn(curIdx, :) = curSpecClassConn;
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
    curClassConn = ...
        conn.axonMeta.axonClass(curSynapses.axonId);
    curClassConn = accumarray( ...
        double(curClassConn), 1, [numel(synTypes), 1]);
    
    curSpecClassConn = ...
        specConn.axonMeta.axonClass(curSynapses.axonId);
    curSpecClassConn = accumarray( ...
        double(curSpecClassConn), 1, [numel(specSynTypes), 1]);
    
    wcT.synapses{curIdx} = curSynapses;
    wcT.classConn(curIdx, :) = curClassConn;
    wcT.specClassConn(curIdx, :) = curSpecClassConn;
end

%% Render isosurface
%{
% Since Mr. Amira is on vacation.
clear cur*;

% Configuration
curCellId = 21;
curCamPos = [44683, 16045, -15604];
curCamTarget = [2852.5, 4320.0, 1965.4];
curCamUpVector = [-0.0178, -0.0785, -0.0153];
curCamViewAngle = 10.4142;

curIso = load(fullfile(isoDir, sprintf('iso-%d.mat', curCellId)));
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
%}

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
dendT.specClassConn = nan(height(dendT), numel(specSynTypes));
dendT.tcExcRatio = nan(size(dendT.id));
dendT.wcRelTcExcRatio = nan(size(dendT.id));

for curIdx = 1:size(dendT, 1)
    curSynapses = dendT.synapses{curIdx};
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
    
    curClassConn = ...
        conn.axonMeta.axonClass(curSynapses.axonId);
    curClassConn = accumarray( ...
        double(curClassConn), 1, [numel(synTypes), 1]);
    
    curSpecClassConn = ...
        specConn.axonMeta.axonClass(curSynapses.axonId);
    curSpecClassConn = accumarray( ...
        double(curSpecClassConn), 1, [numel(specSynTypes), 1]);
    
    dendT.dir(curIdx, :) = curDendDir;
    dendT.classConn(curIdx, :) = curClassConn;
    dendT.specClassConn(curIdx, :) = curSpecClassConn;
    
    dendT.tcExcRatio(curIdx) = ...
        curClassConn(2) / sum(curClassConn(1:2));
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
    fitlm(curDirs, curFracs)
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

%% Correlate input ratios
clear cur*;

curWcMinSyn = 100;
curDendMinSyn = 50;
curRenormNames = {'Excitatory', 'Inhibitory'};

curLims = [0, 1];
curTicks = 0:0.1:1;

curWcT = wcT(sum(wcT.classConn, 2) >= curWcMinSyn, :);
curDendT = dendT(sum(dendT.classConn, 2) >= curDendMinSyn, :);

curData = {curWcT, curDendT};
curSummaries = cell(size(curData));
curSummaryVars = { ...
    'name', ...
    'fit', 'slope', 'relSlope', 'pValue', ...
    'nullModel', 'nullCorrFit', 'nullCorrPValue'};

curDataNames = {'neurons', 'dendrites'};
curVarNames = strcat('frac', categories(specSynTypes));
curVarNames = cat(1, {'ieRatio'}, curVarNames);

for curDataIdx = 1:numel(curData)
    curDataName = curDataNames{curDataIdx};
    curDataT = curData{curDataIdx};
    
    curDataT.tcRatio = ...
        curDataT.classConn(:, 2) ...
     ./ sum(curDataT.classConn(:, 1:2), 2);

    curDataT.ieRatio = ...
        curDataT.classConn(:, 3) ...
     ./ sum(curDataT.classConn(:, 1:3), 2);
    
    for curTypeIdx = 1:numel(specSynTypes)
        curTypeName = char(specSynTypes(curTypeIdx));
        curVarName = sprintf('frac%s', curTypeName);
        
        % NOTE(amotta): Exploit that inhibitory specificity classes are
        % completely independent of the TC / (TC + CC), and that we can
        % derive the expected composition of excitatory specificity classes
        % under the null hypothesis.
        %   In preparation for that, let's renormalize the compositions
        % within the excitatory and inhibitory synapse populations,
        % respectively.
        curRenormIdx = find(cellfun( ...
            @(p) startsWith(curTypeName, p), curRenormNames));
        
        if isempty(curRenormIdx)
            curRenorm = sum(curDataT.specClassConn, 2);
        else
            assert(isscalar(curRenormIdx));
            curRenorm = curRenormNames{curRenormIdx};
            curRenorm = startsWith(categories(specSynTypes), curRenorm);
            curRenorm = sum(curDataT.specClassConn(:, curRenorm), 2);
        end
        
        curDataT.(curVarName) = ...
            curDataT.specClassConn(:, curTypeIdx) ./ curRenorm;
    end
    
    curSummaryT = cell( ...
        numel(curVarNames), ...
        numel(curSummaryVars));
    
    for curVarIdx = 1:numel(curVarNames)
        curVarName = curVarNames{curVarIdx};
        
        curX = curDataT.tcRatio;
        curY = curDataT.(curVarName);
        
        switch curVarName
            case 'ieRatio'
                curNullModel = @(x) ones(size(x)) ...
                 .* synTypeFracs(synTypes == 'Inhibitory') ...
                 ./ sum(synTypeFracs(synTypes ~= 'Other'));
             
            case 'fracOther'
                curNullModel = @(x) ones(size(x)) ...
                 .* synTypeFracs(synTypes == 'Other');
             
            case curVarNames(startsWith(curVarNames, 'fracExcitatory'))
                curExc = 'Excitatory';
                curCcTc = {'Corticocortical', 'Thalamocortical'};
                
                curVar = specSynTypes == curVarName(5:end);
                curRenorm = startsWith(categories(specSynTypes), curExc);
               [~, curIds] = ismember(curCcTc, synTypes);
                curFracs = specAxonTypeFracs(curIds, :);
                
                curNullModel = @(x) feval( ...
                    @(m) ...
                        sum(m(:, curVar, :), 1) ...
                     ./ sum(sum(m(:, curRenorm, :), 1), 2), ...
                    curFracs .* ([1; 0] + [-1; +1] .* x));
                curNullModel = @(x) reshape( ...
                    curNullModel(reshape(x, 1, 1, [])), size(x));
                
            case curVarNames(startsWith(curVarNames, 'fracInhibitory'))
                curInh = 'Inhibitory';
                curVar = specSynTypes == curVarName(5:end);
                curRenorm = startsWith(categories(specSynTypes), curInh);
                
                curFrac = specAxonTypeFracs(synTypes == curInh, :);
                curFrac = curFrac(curVar) / sum(curFrac(curRenorm));
                
                curNullModel = @(x) ones(size(x)) .* curFrac;
                
            otherwise
                error('No null model for variable "%s"', curVarName);
        end
        
        curFit = fitlm(curX, curY);
        curSlope = curFit.Coefficients.Estimate(2);
        curRelSlope = curSlope / mean(curY);
        curPvalue = curFit.Coefficients.pValue(2);
        
        % p value against our custom null model
        curNullCorrY = curY - curNullModel(curX);
        curNullCorrFit = fitlm(curX, curNullCorrY);
        curNullCorrPvalue = curNullCorrFit.Coefficients.pValue(2);
        
        curSummaryT(curVarIdx, :) = { ...
            curVarName, ...
            curFit, curSlope, curRelSlope, curPvalue, ...
            curNullModel, curNullCorrFit, curNullCorrPvalue};
    end
    
    curSummaryT = cell2table( ...
        curSummaryT, 'VariableNames', curSummaryVars);
    
    fprintf('Results for %s\n', curDataName);
    curShowT = {'fit', 'nullModel', 'nullCorrFit'};
    curShowT = setdiff(curSummaryVars, curShowT, 'stable');
    disp(sortrows(curSummaryT(:, curShowT), 'relSlope', 'ascend'));
    
    curData{curDataIdx} = curDataT;
    curSummaries{curDataIdx} = curSummaryT;
end

% Plot correlation of dendrite length with TC input fraction
curDataT = strcmpi(curDataNames, 'dendrites');
curDataT = curData{curDataT};

curX = cellfun(@max, curDataT.nodeDists) / 1E3;
curY = curDataT.tcRatio;

curDistLims = [0, ceil(max(curX) / 10) * 10];
curDistTicks = curDistLims(1):50:curDistLims(end);

curFit = fitlm(curX, curY);
curFitY = curFit.predict(curDistLims(:))';

curFig = figure();
curAx = axes(curFig);
hold(curAx, 'on');

curScatter = scatter(curAx, curX, curY, 20, 'o');
curScatter.MarkerEdgeColor = 'none';
curScatter.MarkerFaceColor = colors(1, :);

curFitPlot = plot( ...
    curAx, curDistLims, curFitY, ...
    'Color', 'black', 'LineWidth', 2);

axis(curAx, 'square');
xlim(curAx, curDistLims);
xticks(curAx, curDistTicks);
ylim(curAx, [0, 1]);

xlabel(curAx, 'Dendrite length [µm]');
ylabel(curAx, 'tcRatio');

curLeg = { ...
    sprintf( ...
        'Dendrites (n = %d)', height(curDataT)), ...
    sprintf( ...
        'y = %.2f %+.2fx (p_{slope} = %.2f)', ...
        curFit.Coefficients.Estimate, ...
        curFit.Coefficients.pValue(end))};
curLeg = legend(curLeg);

connectEM.Figure.config(curFig, info);
curFig.Position(3:4) = [360, 360];

% Plot correlates of TC input fraction
curFig = figure();
for curDataIdx = 1:numel(curData)
    curDataName = curDataNames{curDataIdx};
    curSummaryT = curSummaries{curDataIdx};
    curDataT = curData{curDataIdx};
    curX = curDataT.tcRatio;
    
    for curVarIdx = 1:numel(curVarNames)
        curVarName = curVarNames{curVarIdx};
        curY = curDataT.(curVarName);
        
        curFit = curSummaryT.fit{curVarIdx};
        curFitY = curFit.predict(curLims(:))';
        
        curNull = curSummaryT.nullModel{curVarIdx};
        curNullY = curNull(curLims);
        
        curNullCorrFit = curSummaryT.nullCorrFit{curVarIdx};
        
        curAx = numel(curVarNames) * (curDataIdx - 1) + curVarIdx;
        curAx = subplot(numel(curData), numel(curVarNames), curAx);
        hold(curAx, 'on');

        curScatter = scatter(curAx, curX, curY, 20, 'o');
        curScatter.MarkerEdgeColor = 'none';
        curScatter.MarkerFaceColor = colors(1, :);

        curFitPlot = plot( ...
            curAx, curLims, curFitY, ...
            'Color', 'black', 'LineWidth', 2);
        curNullPlot = plot( ...
            curAx, curLims, curNullY, ...
            'Color', 'black', 'LineStyle', '--');

        axis(curAx, 'equal');
        xlim(curAx, curLims); xticks(curAx, curTicks);
        ylim(curAx, curLims); yticks(curAx, curTicks);
        
        xlabel(curAx, sprintf('%s (%s)', 'tcRatio', curDataName));
        ylabel(curAx, sprintf('%s (%s)', curVarName, curDataName));
        
        curLeg = { ...
            sprintf( ...
                'y = %.2f %+.2fx (p_{slope} = %.2f)', ...
                curFit.Coefficients.Estimate, ...
                curFit.Coefficients.pValue(end)), ...
            sprintf( ...
                'Null model (p_{slope vs. null} = %.2f)', ...
                curNullCorrFit.Coefficients.pValue(end))};
        curLeg = legend([curFitPlot, curNullPlot], curLeg);
        curLeg.Location = 'north';
    end
end

connectEM.Figure.config(curFig, info);
curFig.Position(3:4) = [6600, 1000];

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
        curFig = figure();
    elseif ~isempty(plotDir) %#ok
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
    
    curCellCount = numel(unique(curSyns.parentId)) %#ok
    curDendCount = sum(ismember(dendT.id, curSyns.parentId)) %#ok
    curSynCount = height(curSyns) %#ok
    
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
