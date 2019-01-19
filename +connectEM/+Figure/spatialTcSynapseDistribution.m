% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

maxSamples = 100;
bandWidth = 5;
marginUm = 3;

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

%% Calculate TC synapse density volume
clear cur;
gridSize = max(limits(:, 2));
gridSize = round(maxSamples .* limits(:, 2) ./ gridSize);

grid = cell(1, 1, 1, 3);
[grid{:}] = ndgrid( ...
    linspace(limits(1, 1), limits(1, 2), gridSize(1)), ...
    linspace(limits(2, 1), limits(2, 2), gridSize(2)), ...
    linspace(limits(3, 1), limits(3, 2), gridSize(3)));
grid = reshape(cell2mat(grid), [], 3);

curSynT = synT;
curSynT = curSynT( ...
    all(curSynT.pos > limits(:, 1)', 2) ...
  & all(curSynT.pos < limits(:, 2)', 2), :);

tic;
densVol = mvksdensity( ...
    curSynT.pos, grid, ...
    'Support', transpose(limits), ...
    'BoundaryCorrection', 'reflection', ...
    'BandWidth', bandWidth);
toc;
densVol = reshape(densVol, gridSize(:)');

%% Try isosurface visualization
fig = figure;
ax = axes(fig);
hold(ax, 'on');

curColors = jet();

quantVec = [0.5, 0.85];
for curQuant = quantVec
    curThresh = quantile(densVol(:), curQuant);
    curIso = isosurface(densVol, curThresh);
    curIso.vertices = curIso.vertices(:, [2, 1, 3]);
    curPatch = patch(curIso);
    curPatch.EdgeColor = 'none';
    
    curColor = curQuant / max(quantVec);
    curColor = ceil(curColor * size(curColors, 1));
    curPatch.FaceColor = curColors(curColor, :);
    curPatch.FaceAlpha = curQuant / 2;
end

axis(ax, 'equal');
xlim(ax, [1, size(densVol, 1)]);
ylim(ax, [1, size(densVol, 2)]);
zlim(ax, [1, size(densVol, 3)]);
set(ax, 'Box', 'on', 'BoxStyle', 'full');

view(ax, 3);
camlight(ax);

%% Tangential view
clear cur*;
fig = figure();
fig.Color = 'white';

curImData = shiftdim(sum(densVol, 1), 1);
curImData = uint8(double(intmax('uint8')) * curImData / max(curImData(:)));
curImRatio = size(curImData, 1) / size(curImData, 2);

ax = axes(fig);
curIm = image(ax, curImData);
colormap(ax, jet(256));
axis(ax, 'image');

ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
xlim(ax, [1, size(curImData, 2)]);
ylim(ax, [1, size(curImData, 1)]);

title(ax, { ...
    info.filename; info.git_repos{1}.hash; ...
    'Tangential projections of TC synapse density'}, ...
    'FontWeight', 'normal', 'FontSize', 10);
