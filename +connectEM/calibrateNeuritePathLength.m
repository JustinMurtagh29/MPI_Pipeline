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
trunkLUT = Agglo.buildLUT(maxSegId, trunks);

%% Loading dendrite calibration data
dendCalibT = ...
    connectEM.Dendrite.Data.getDir('pathLengthCalibration');
dendCalibT = ...
    connectEM.Connectome.loadPathLengthCalibrationNmls(param, dendCalibT);

dendCalibT.trunkId = cellfun( ...
    @(segIds) mode(nonzeros(trunkLUT(segIds))), ...
    conn.dendrites(dendCalibT.id));
dendCalibT.trunkPathLength = ...
    lengths.trunkPathLengths(dendCalibT.trunkId);

%% Calculating correction coefficient
corrCoeff = mean(dendCalibT.pathLength ./ dendCalibT.trunkPathLength);

%% Plot
curFig = figure();
curFig.Color = 'white';
curFig.Position(3:4) = [420, 420];

curAx = axes(curFig);
axis(curAx, 'equal');
hold(curAx, 'on');

scatter(curAx, ...
    dendCalibT.trunkPathLength / 1E3, ...
    dendCalibT.pathLength / 1E3, 72, '.');

curLim = max(curAx.XLim(2), curAx.YLim(2));
curAx.XLim = [0, curLim];
curAx.YLim = [0, curLim];

plot(curAx, curAx.XLim, corrCoeff * curAx.XLim);
plot(curAx, curAx.XLim, curAx.YLim, 'Color', 'black', 'LineStyle', '--');

curAx.TickDir = 'out';
xlabel(curAx, 'MST-based trunk length (µm)');
ylabel(curAx, 'Tracing-based trunk length (µm)');

curLeg = legend(curAx, { ...
    'Calibration points', ...
    sprintf('Correction (coefficient = %.3f)', corrCoeff)}, ...
    'Location', 'NorthWest');
curLeg.Box = 'off';

title(curAx,  ...
    {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
