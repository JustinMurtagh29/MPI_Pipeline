% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');
trunkFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3.mat');
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto.mat');

[~, lengthFile] = fileparts(connFile);
lengthFile = sprintf('%s_pathLengths.mat', lengthFile);
lengthFile = fullfile(fileparts(connFile), lengthFile);

calibClasses = { ...
    'AxonInitialSegment', ...
    'ApicalDendrite', ...
    'SmoothDendrite', ...
    'Random'};

plotTargetClasses = { ...
    'AxonInitialSegment', 'AIS'; ...
    'ApicalDendrite', 'AD'; ...
    'SmoothDendrite', 'SD'; ...
    'OtherDendrite', 'Other'};

plotLabels = plotTargetClasses(:, 2);
plotTargetClasses = plotTargetClasses(:, 1);

minSynPost = 10;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

conn = connectEM.Connectome.load(param, connFile);

trunks = load(trunkFile);
trunks = trunks.dendrites(trunks.indBigDends);
trunks = Agglo.fromSuperAgglo(trunks);

shAgglos = load(shFile, 'shAgglos');
shAgglos = shAgglos.shAgglos;

lengths = load(lengthFile);

maxSegId = Seg.Global.getMaxSegId(param);

%% Find trunk length for dendrite in connectome
trunkLUT = Agglo.buildLUT(maxSegId, trunks);
conn.denMeta.trunkId = cellfun( ...
    @(ids) mode(nonzeros(trunkLUT(ids))), conn.dendrites);

curMask = ~isnan(conn.denMeta.trunkId);
conn.denMeta.trunkLength = nan(height(conn.denMeta), 1);
conn.denMeta.trunkLength(curMask) = ...
    lengths.trunkPathLengths(conn.denMeta.trunkId(curMask));
clear curMask;

%%
dendT = conn.denMeta;
dendT = dendT( ...
    dendT.targetClass == 'AxonInitialSegment' ...
  | dendT.synCount >= minSynPost, :);
dendT = dendT(ismember(dendT.targetClass, plotTargetClasses), :);

dendT.spineCount = ...
    connectEM.Dendrite.calculateSpineCount( ...
        param, conn.dendrites(dendT.id), shAgglos);

%% Process NML files
calibData = struct;
for curCalibClassIdx = 1:numel(calibClasses)
    curCalibClass = calibClasses{curCalibClassIdx};
    
    curDir = sprintf('pathLengthCalibration%s', curCalibClass);
    curDir = connectEM.Dendrite.Data.getFile(curDir);
    
    curNmlFiles = dir(fullfile(curDir, '*.nml'));
    curNmlFiles = reshape({curNmlFiles.name}, [], 1);
    curNmlFiles = fullfile(curDir, curNmlFiles);

    curTaskFile = fullfile(curDir, 'taskIds.txt');
    curUseTaskIds = true;
    
    try
        curTasks = connectEM.Chiasma.Util.loadTaskIds(curTaskFile);
    catch
        curTasks = table;
        curTasks.nmlFile = curNmlFiles;
        curUseTaskIds = false;
    end
    
    % Extract agglomerate ID
    curTasks.aggloId = regexpi( ...
        curTasks.nmlFile, '.*-(\d+)\.nml', 'tokens', 'once');
    assert(all(cellfun(@isscalar, curTasks.aggloId)));
    
    curTasks.aggloId = cat(1, curTasks.aggloId{:});
    curTasks.aggloId = cellfun(@str2double, curTasks.aggloId);
    
    if curUseTaskIds
        second = @(in) in(2);
        curNmlTaskId = cellfun(@(name) ...
            second(strsplit(name, '__')), curNmlFiles);

       [~, curTasks.nmlFile] = ismember(curTasks.id, curNmlTaskId);
        
        curTasks(~curTasks.nmlFile, :) = [];
        curTasks.nmlFile = curNmlFiles(curTasks.nmlFile);
    end

    % Calculate tracing-based path length
    curTasks.calibLength = nan(height(curTasks), 1);

    for curTaskIdx = 1:height(curTasks)
        curNmlFile = curTasks.nmlFile{curTaskIdx};
        curSkel = skeleton(curNmlFile);

        % Sanity check
        assert(curSkel.numTrees() == 2);
        assert(any(curSkel.thingIDs == 1));

        curCalibTreeId = find(curSkel.thingIDs ~= 1);
        curCalibLength = curSkel.pathLength( ...
            curCalibTreeId, param.raw.voxelSize);

        curTasks.calibLength(curTaskIdx) = curCalibLength;
    end
    
    calibData.(curCalibClass) = curTasks;
end

%% Visualize results
clear cur*;

curFig = figure();
curFig.Color = 'white';
curFig.Position(3:4) = [1025, 350];

for curIdx = 1:numel(calibClasses)
    curCalibClass = calibClasses{curIdx};
    curCalibT = calibData.(curCalibClass);
    
    curCalibT.autoLength = conn.denMeta.trunkLength(curCalibT.aggloId);
    
    % NOTE(amotta): This is what we did previously. This fits both an
    % offset and a slope parameter. However, for the axon and dendrite path
    % length calibration (connectEM.calibrateNeuritePathLength). we did not
    % allow non-zero offsets. So, let's use the same approach here.
    %
    % curInterp = fit(curCalibT.autoLength, curCalibT.calibLength, 'poly1');
    
    curInterp = fit( ...
        [0; sum(curCalibT.autoLength)], ...
        [0; sum(curCalibT.calibLength)], 'poly1');

    curAx = subplot(1, numel(calibClasses), curIdx);
    axis(curAx, 'square');
    hold(curAx, 'on');

    scatter(curAx, ...
        curCalibT.autoLength / 1E3, ...
        curCalibT.calibLength / 1E3, ...
        80, '.');

    curLimits = [0, max(curAx.XLim(2), curAx.YLim(2))];
    curAx.XLim = curLimits; curAx.YLim = curLimits;

    plot(curAx, curLimits(:), curInterp(1E3 * curLimits) / 1E3);
    plot(curAx, curLimits, curLimits, 'Color', 'black', 'LineStyle', '--');
    
    curTitle = { ...
        curCalibClass; ...
        sprintf('%d calibration points', height(curCalibT)); ...
        sprintf('Correction coefficient: %.3f', curInterp.p1)};
    title(curAx, curTitle, 'FontWeight', 'normal', 'FontSize', 10);
end

curAxes = flip(curFig.Children);
set(curAxes, 'TickDir', 'out');

xlabel(curAxes(1), 'MST-based trunk length (µm)');
ylabel(curAxes(1), 'Tracing-based trunk length (µm)');

annotation( ...
    curFig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
	'String', {info.filename; info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');

%% Correct path length
dendT.corrTrunkLength = nan(size(dendT.trunkLength));
for curIdx = 1:numel(plotTargetClasses)
    curTargetClass = plotTargetClasses{curIdx};
    curCalibClass = calibClasses{curIdx};
    
    curCalibT = calibData.(curCalibClass);
    curCalibT.autoLength = conn.denMeta.trunkLength(curCalibT.aggloId);
    
    % NOTE(amotta): This is what we used to do here. See above for details.
    % curInterp = fit(curCalibT.autoLength, curCalibT.calibLength, 'poly1');
    
    curInterp = fit( ...
        [0; sum(curCalibT.autoLength)], ...
        [0; sum(curCalibT.calibLength)], 'poly1');
    
    curMask = dendT.targetClass == curTargetClass;
    dendT.corrTrunkLength(curMask) = curInterp(dendT.trunkLength(curMask));
end

dendT.correSpineDensity = ...
    dendT.spineCount ./ ( ...
    dendT.corrTrunkLength / 1E3);

%% Plot
pathEdges = linspace(0, 200, 41);
pathClasses = [1, 3, 2];

densityEdges = linspace(0, 2.5, 51);
densityClasses = [1, 3, 2, 4];

% Make sure that `pathClasses` and `densityClasses` share prefix.
assert(feval( ...
    @(n) isequal(pathClasses(1:n), densityClasses(1:n)), ...
    min(numel(pathClasses), numel(densityClasses))));

colors = get(groot, 'defaultAxesColorOrder');

plotHist = @(ax, edges, data) ...
    histogram( ...
        ax, data, ...
        'BinEdges', edges, ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2, ...
        'FaceAlpha', 1);

fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [750, 365];

% Path length
ax = subplot(1, 2, 1);
hold(ax, 'on');

for curTargetClassId = flip(pathClasses)
    curTargetClass = plotTargetClasses(curTargetClassId);
    curColor = colors(curTargetClassId, :);
    
    curData = dendT.targetClass == curTargetClass;
    curData = dendT.corrTrunkLength(curData) / 1E3;
    curHist = plotHist(ax, pathEdges, curData);
    curHist.EdgeColor = curColor;
end

xlim(ax, pathEdges([1, end]));
xlabel(ax, 'Path length (µm)');
ylabel(ax, 'Dendrites');

% Spine density
ax = subplot(1, 2, 2);
hold(ax, 'on');

for curTargetClassId = flip(densityClasses)
    curTargetClass = plotTargetClasses(curTargetClassId);
    curColor = colors(curTargetClassId, :);
    
    curData = dendT.targetClass == curTargetClass;
    curData = dendT.correSpineDensity(curData, :);
    curHist = plotHist(ax, densityEdges, curData);
    curHist.EdgeColor = curColor;
end

xlim(ax, densityEdges([1, end]));
xlabel(ax, 'Spine density (µm^{-1})');

set(fig.Children, ...
    'TickDir', 'out', ...
    'PlotBoxAspectRatio', ones(1, 3), ...
    'DataAspectRatioMode', 'auto');

% NOTE(amotta): `pathClasses` and `densityClasses` might be different. But
% they share the same prefix (see above sanity check). So, we can just use
% the longer of the two vectors. This is equivalent to the `stable` union.
leg = flip(union(pathClasses, densityClasses, 'stable'));
leg = legend(ax, plotLabels(leg), 'Location', 'NorthEast');
leg.Box = 'off';

annotation( ...
    fig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
	'String', {info.filename; info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');
