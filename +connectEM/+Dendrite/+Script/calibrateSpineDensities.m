% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');
trunkFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3.mat');

[~, lengthFile] = fileparts(connFile);
lengthFile = sprintf('%s_pathLengths.mat', lengthFile);
lengthFile = fullfile(fileparts(connFile), lengthFile);

targetClasses = { ...
    'ApicalDendrite', ...
    'SmoothDendrite', ...
    'AxonInitialSegment'};

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

conn = connectEM.Connectome.load(param, connFile);

trunks = load(trunkFile);
trunks = trunks.dendrites(trunks.indBigDends);
trunks = Agglo.fromSuperAgglo(trunks);

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

%% Process NML files
calibData = cell(size(targetClasses));
for curTargetClassIdx = 1:numel(targetClasses)
    curTargetClass = targetClasses{curTargetClassIdx};
    
    curDir = sprintf('pathLengthCalibration%s', curTargetClass);
    curDir = connectEM.Dendrite.Data.getFile(curDir);

    curTaskFile = fullfile(curDir, 'taskIds.txt');
    curTasks = connectEM.Chiasma.Util.loadTaskIds(curTaskFile);

    % Extract agglomerate ID
    curTasks.aggloId = regexpi( ...
        curTasks.nmlFile, '.*-(\d+)\.nml', 'tokens', 'once');
    assert(all(cellfun(@isscalar, curTasks.aggloId)));

    curTasks.aggloId = cat(1, curTasks.aggloId{:});
    curTasks.aggloId = cellfun(@str2double, curTasks.aggloId);

    % Find NML file with results
    curNmlFiles = dir(fullfile(curDir, '*.nml'));
    curNmlFiles = reshape({curNmlFiles.name}, [], 1);

    second = @(in) in(2);
    curNmlTaskId = cellfun(@(name) ...
        second(strsplit(name, '__')), curNmlFiles);

    [~, curTasks.nmlFile] = ismember(curTasks.id, curNmlTaskId);
    curTasks(~curTasks.nmlFile, :) = [];

    curTasks.nmlFile = curNmlFiles(curTasks.nmlFile);
    curTasks.nmlFile = fullfile(curDir, curTasks.nmlFile);

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
    
    calibData{curTargetClassIdx} = curTasks;
end

%% Visualize results
clear cur*;

curFig = figure();
curFig.Color = 'white';
curFig.Position(3:4) = [1025, 350];

for curIdx = 1:numel(targetClasses)
    curTargetClass = targetClasses{curIdx};
    
    curCalibT = calibData{curIdx};
    curCalibT.autoLength = conn.denMeta.trunkLength(curCalibT.aggloId);
    curInterp = fit(curCalibT.autoLength, curCalibT.calibLength, 'poly1');

    curAx = subplot(1, numel(targetClasses), curIdx);
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
    
    title(curAx, curTargetClass, 'FontWeight', 'normal', 'FontSize', 10);
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
