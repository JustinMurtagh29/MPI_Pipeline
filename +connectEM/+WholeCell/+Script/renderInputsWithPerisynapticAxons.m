% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');
wcIsoDir = '/tmpscratch/amotta/l4/2018-05-10-whole-cell-isosurfaces-spine-evolution/full';
periSynIsoDir = '/tmpscratch/amotta/l4/2019-10-20-whole-cell-input-axon-snippet-isosurfaces_v2';

colors = get(groot, 'defaultAxesColorOrder');
colors = colors(1:4, :);

colorClasses = { ...
    'Corticocortical', 'Thalamocortical', ...
    'Inhibitory', 'Other'};
assert(isequal(numel(colorClasses), size(colors, 1)));

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, syn, axonClasses] = connectEM.Connectome.load(param, connFile);
synT = connectEM.Connectome.buildSynapseTable(conn, syn);

periSynT = fullfile(periSynIsoDir, 'input-synapses.mat');
periSynT = Util.load(periSynT, 'synT');

%% Render cells
clear cur*;

for curCellId = 5 % [1, 5, 21]
    curDendIds = find(conn.denMeta.cellId == curCellId);
    curSynT = synT(ismember(synT.postAggloId, curDendIds), :);
    
    curSynT.colorId = conn.axonMeta.axonClass(curSynT.preAggloId);
   [~, curSynT.colorId] = ismember(curSynT.colorId, colorClasses);
    
   [~, curSynT.isoId] = ismember(curSynT.id, periSynT.id);
    curSynT = curSynT(logical(curSynT.isoId), :);
    
    curFig = figure();
    curAx = axes(curFig); %#ok
    
    hold(curAx, 'on');
    daspect(curAx, 1 ./ param.raw.voxelSize);
    
    % Neuron
    curCellIso = loadIsoSurf(wcIsoDir, curCellId);
    curCellIso = patch(curAx, curCellIso);
    curCellIso.EdgeColor = 'none';
    curCellIso.FaceColor = repelem(0.85, 3);
    material(curCellIso, 'dull');
    
    tic;
    rng(0);
    for curSynIdx = 1:height(curSynT)
        curSynIsoId = curSynT.isoId(curSynIdx);
        curSynColor = colors(curSynT.colorId(curSynIdx), :);
        
        % Randomize color
        curSynColor = rgb2hsv(curSynColor);
        curSynColor(1) = curSynColor(1) + 0.1 * (rand() - 0.5);
        curSynColor(1) = mod(curSynColor(1), 1);
        curSynColor = hsv2rgb(curSynColor);
        
        curSynIso = loadIsoSurf(periSynIsoDir, curSynIsoId);
        curSynIso = patch(curAx, curSynIso);
        curSynIso.FaceColor = curSynColor;
        curSynIso.EdgeColor = 'none';
        material(curSynIso, 'dull');
        
        Util.progressBar(curSynIdx, height(curSynT));
    end
    
    curFig.Color = 'none';
    curAx.Visible = 'off';
    curAx.Color = 'none';
    
    view(curAx, -90, 90);
    camlight(curAx);
end

%% Utilities
function isoSurf = loadIsoSurf(dir, id)
    isoSurf = fullfile(dir, 'mat', sprintf('iso-%d.mat', id));
    isoSurf = Util.load(isoSurf, 'isoSurf');
end
