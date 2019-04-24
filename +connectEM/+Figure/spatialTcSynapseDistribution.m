% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

imSize = 128;
bandWidth = 10;
marginUm = 3;
dimIds = [2, 3];

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, syn] = connectEM.Connectome.load(param, connFile);

%% Preparing data
synT = connectEM.Connectome.buildSynapseTable(conn, syn);
synT.synType = conn.axonMeta.axonClass(synT.preAggloId);
synT = synT(synT.synType == 'Thalamocortical', :);

synPos = connectEM.Synapse.calculatePositions(param, syn);
synT.pos = synPos(synT.id, :);
clear synPos;

% Physical units, relative to center
synT.pos = synT.pos - mean(param.bbox, 2)';
synT.pos = synT.pos .* param.raw.voxelSize ./ 1E3;

%% Determine bounding box
clear cur*;
limits = param.bbox(:, 2) - mean(param.bbox, 2);
limits = limits .* param.raw.voxelSize(:) / 1E3;
limits = limits - marginUm;
limits = [-1, +1] .* limits;

%% Tangential view
clear cur*;

curSynT = synT;
curSynT = curSynT( ...
    all(curSynT.pos > limits(:, 1)', 2) ...
  & all(curSynT.pos < limits(:, 2)', 2), :);

[~, curImData] = connectEM.Libs.kde2d( ...
    curSynT.pos(:, dimIds), imSize, ...
    limits(dimIds, 1)', limits(dimIds, 2)', ...
    repelem(bandWidth, 1, 2));

% NOTE(amotta): Let's convert the density estimate into an estimate of the
% spatial synapse density. I assume that there exists an elegant analytical
% way to calculate this. Instead, let's just make sure that the mean value
% is right (and zero is zero).
curImData = ...
    curImData / mean(curImData(:)) ...
  * size(curSynT, 1) / prod(diff(limits, 1, 2));

curImData = transpose(curImData);
curImMax = max(curImData(:));

fig = figure();
fig.Color = 'white';

ax = axes(fig);
hold(ax, 'on');

curIm = imagesc(ax, curImData);
ax.CLim = [0, curImMax];

curIm.XData = limits(dimIds(1), :);
curIm.YData = limits(dimIds(2), :);

colormap(ax, jet(256));
colorbar('peer', ax);

curCont = linspace(0, curImMax, 5);
curCont = curCont(2:(end - 1));

[~, curCont] = contour(ax, ...
    double(curImData), curCont, ...
    'LineColor', 'black');

curCont.XData = ...
    limits(dimIds(1), 1) ...
  + diff(limits(dimIds(1), :)) ...
  * (curCont.XData - 1) / (imSize - 1);
curCont.YData = ...
    limits(dimIds(2), 1) ...
  + diff(limits(dimIds(2), :)) ...
  * (curCont.YData - 1) / (imSize - 1);

axis(ax, 'image');

ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
xlim(ax, limits(dimIds(1), :));
ylim(ax, limits(dimIds(2), :));

connectEM.Figure.config(fig, info);
fig.Position(3:4) = [250, 180];
