% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');
outputMapFile = '/tmpscratch/amotta/l4/2018-07-26-tracing-based-output-maps/20180726T190355_results.mat';
isoDir = '/tmpscratch/amotta/l4/2018-05-10-whole-cell-isosurfaces-spine-evolution/full/mat/';

% From an email by HW on 04.06.2018
cellPos = [236, 3498, 990];

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

segCentroid = Seg.Global.getSegToCentroidMap(param);
segMass = Seg.Global.getSegToSizeMap(param);

[conn, syn] = connectEM.Connectome.load(param, connFile);
syn.synapses.pos = connectEM.Synapse.calculatePositions(param, syn);

%% Find cell ID
somata = conn.denMeta;
somata(somata.targetClass ~= 'Somata', :) = [];

wmean = @(w, v) sum((w / sum(w)) .* v, 1);
somata.pos = cell2mat(cellfun( ...
    @(ids) wmean(segMass(ids), segCentroid(ids, :)), ...
    conn.dendrites(somata.id), 'UniformOutput', false));

[~, cellId] = pdist2(somata.pos, cellPos, 'euclidean', 'Smallest', 1);
cellId = somata.cellId(cellId);

%% Load data for visualization
iso = load(fullfile(isoDir, sprintf('iso-%d.mat', cellId)));
iso = iso.isoSurf;

outputMap = load(outputMapFile);

%% Find corresponding NML file
clear cur*;

curNmlFiles = {outputMap.axonData.nmlFile};
curNmlSomaPos = nan(numel(curNmlFiles), 3);

for curIdx = 1:numel(curNmlFiles)
    curNml = skeleton(curNmlFiles{curIdx});
    curNmlDendId = find(strcmpi(curNml.names, 'Dendrite'));
    curNmlSomaId = curNml.getNodesWithComment('Soma', curNmlDendId, 'insensitive');
    curNmlSomaPos(curIdx, :) = curNml.nodes{curNmlDendId}(curNmlSomaId, 1:3);
end

[~, outputMapId] = pdist2( ...
    curNmlSomaPos, cellPos, ...
    'euclidean', 'Smallest', 1);

synIds = outputMap.axonData(outputMapId).synapses.id;
nmlFile = curNmlFiles{outputMapId};

%% Load data for visualization
nml = slurpNml(nmlFile);
tree = NML.buildTreeTable(nml);
tree = tree(strcmpi(tree.name, 'Axon'), :);

axon = struct;
axon.nodes = [ ...
    tree.nodes{1}.x, ...
    tree.nodes{1}.y, ...
    tree.nodes{1}.z] + 1;
axon.edges = [ ...
    tree.edges{1}.source, ...
    tree.edges{1}.target];
[~, axon.edges] = ismember( ...
    axon.edges, tree.nodes{1}.id);

%% Generate isosurface
fig = figure();
fig.Color = 'none';

ax = axes(fig);
hold(ax, 'on');
view(ax, 90, 90);

ax.Color = 'none';
ax.Visible = 'off';

daspect(ax, 1 ./ param.raw.voxelSize);

p = patch(ax, iso);
p.EdgeColor = 'none';
p.FaceColor = 'blue';
material(p, 'dull');
camlight(ax);

% NOTE(amotta): Here, I'm concatenating the X and Y values of all
% edges, respectively. To break up the chain, we have to intersperse
% NaNs between each pair of X or Y values. This allows us to use the
% `plot` command despite the coordinate sequence to be disconnected.
% This is important because calling and displaying the results of the
% `plot` command is dramatically faster than creating `lines` for each
% edge.
curEdges = transpose(axon.edges);
curX = reshape(axon.nodes(curEdges(:), 1), 2, []);
curX = reshape(cat(1, curX, nan(1, size(curX, 2))), 1, []);
curY = reshape(axon.nodes(curEdges(:), 2), 2, []);
curY = reshape(cat(1, curY, nan(1, size(curY, 2))), 1, []);
curZ = reshape(axon.nodes(curEdges(:), 3), 2, []);
curZ = reshape(cat(1, curZ, nan(1, size(curZ, 2))), 1, []);

curPlot = plot3(ax, curX, curY, curZ);
curPlot.Color = 'black';
curPlot.LineWidth = 2;

[curX, curY, curZ] = sphere();
curX = curX / param.raw.voxelSize(1);
curY = curY / param.raw.voxelSize(2);
curZ = curZ / param.raw.voxelSize(3);
curRad = 1.5E3;

synColors = [0, 0, 0; 1, 0, 1];
synIsSpine = syn.synapses.type(synIds) == 'PrimarySpine';
synColors = synColors(1 + synIsSpine, :);
synPos = syn.synapses.pos(synIds, :);

for curId = 1:numel(synIds)
    curColor = synColors(curId, :);
    curPos = synPos(curId, :);
    
    curSurf = surf(ax, ...
        curRad * curX + curPos(1), ...
        curRad * curY + curPos(2), ...
        curRad * curZ + curPos(3));
    curSurf.FaceColor = curColor;
    curSurf.EdgeColor = 'none';
    material(curSurf, 'dull');
end

view(ax, 90, 90);

annotation( ...
    fig, 'textbox', [0, 0.9, 1, 0.1], ...
    'String', {info.filename; info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');
