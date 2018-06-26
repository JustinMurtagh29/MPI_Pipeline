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

%% Find trunk to dendrites
trunkLUT = Agglo.buildLUT(maxSegId, trunks);
conn.denMeta.trunkId = cellfun( ...
    @(segIds) mode(nonzeros(trunkLUT(segIds))), conn.dendrites);

%% Loading axon calibration data
axonCalibT = connectEM.Axon.Data.getDir('pathLengthCalibration');
axonCalibT = connectEM.loadPathLengthCalibrationNmls(param, axonCalibT);

axonCalibT.autoPathLength = lengths.axonPathLengths(axonCalibT.id);

%% Loading dendrite calibration data
dendCalibT = connectEM.Dendrite.Data.getFile('pathLengthCalibrationRandom');
dendCalibT = connectEM.loadPathLengthCalibrationNmls(param, dendCalibT);

dendCalibT.trunkId = conn.denMeta.trunkId(dendCalibT.id);
dendCalibT.autoPathLength = lengths.trunkPathLengths(dendCalibT.trunkId);

%% Calculating correction coefficient
axonCorrCoeff = ...
    sum(axonCalibT.pathLength) ...
  / sum(axonCalibT.autoPathLength);
dendCorrCoeff = ...
    sum(dendCalibT.pathLength) ...
  / sum(dendCalibT.autoPathLength);

%% Plot
curPlotConfigs = { ...
    'Axons', axonCalibT; ...
    'Dendrites', dendCalibT};

curFig = figure();
curFig.Color = 'white';
curFig.Position(3:4) = [800, 440];

for curIdx = 1:size(curPlotConfigs, 1)
    curName = curPlotConfigs{curIdx, 1};
    curCalibT = curPlotConfigs{curIdx, 2};
    
    curCorrCoeff = ...
        sum(curCalibT.pathLength) ...
      / sum(curCalibT.autoPathLength);
    
    curAx = subplot(1, size(curPlotConfigs, 1), curIdx);
    axis(curAx, 'equal');
    hold(curAx, 'on');

    scatter(curAx, ...
        curCalibT.autoPathLength / 1E3, ...
        curCalibT.pathLength / 1E3, 80, '.');

    curLim = max(curAx.XLim(2), curAx.YLim(2));
    curAx.XLim = [0, curLim];
    curAx.YLim = [0, curLim];

    plot(curAx, ...
        curAx.XLim, curCorrCoeff * curAx.XLim);
    plot(curAx, ...
        curAx.XLim, curAx.YLim, ...
        'Color', 'black', 'LineStyle', '--');
    
    title(curAx,{ ...
        curName; ...
        sprintf([ ...
            '%d calibration points. ', ...
            'Correction coefficient: %.3f'], ...
            height(curCalibT), curCorrCoeff)}, ...
        'FontWeight', 'normal', 'FontSize', 10);
end

curAx = curFig.Children(end);
xlabel(curAx, 'MST-based path length (µm)');
ylabel(curAx, 'Tracing-based path length (µm)');

set(curFig.Children, 'TickDir', 'out');

annotation( ...
    curFig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
	'String', {info.filename; info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');

