% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-linearized_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto.mat');

% NOTE(amotta): I've previously inspected a bunch of multi-synaptic
% connections (see notes 2018-11-15-Mono-vs-Bi-Spiny-Connections).
%   Below is a hand-picked subset of these connections that lends itself
% to illustrate the axon-spine interface (e.g., not too oblique, no
% dramatic alignment issues.
%   The name `adsT` stands for `axon dendrite spine-head table`.
adsT = [ ...
    12236,  5946, 258240; % https://webknossos.brain.mpg.de/annotations/Explorational/5bedb18501000028208a2b33
    52083,  3774,  31809; % https://webknossos.brain.mpg.de/annotations/Explorational/5bedb2fd01000035218a2b84
    48188, 11119, 264330; % https://webknossos.brain.mpg.de/annotations/Explorational/5bedb4c8010000f0228a2bdd
    21002,  3873, 339464; % https://webknossos.brain.mpg.de/annotations/Explorational/5bedb88c010000ef238a2cad
    37683,   846,  72632; % https://webknossos.brain.mpg.de/annotations/Explorational/5bedb93a010000ea248a2ccd
];

adsT = array2table( ...
    adsT, 'VariableNames', ...
    {'axonId', 'dendId', 'shId'});

colors = get(groot, 'defaultAxesColorOrder');
colors = colors(1:2, :);

asiColor = colors(2, :);
asiAlpha = 0.9;

% NOTE(amotta): If set to non-NaN values, the user-specified position will
% be used instead of the automatically computed ones.
adsT.pos = nan(height(adsT), 3);
adsT.isoFovUm = 5 * ones(height(adsT), 3);
adsT.emFovUm(:) = 2;

% See https://gitlab.mpcdf.mpg.de/connectomics/amotta/blob/b3d6b7be4f876c54e4261026d5ab47d4edd49d89/matlab/+L4/+Figure/highResEmSamples.m
% Range from Benedikt, used for SynEM paper.
emRange = [60, 180];

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);
points = Seg.Global.getSegToPointMap(param);
sizes = Seg.Global.getSegToSizeMap(param);

[conn, syn, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);

% Loading spine head agglomerates
shAgglos = load(shFile, 'shAgglos');
shAgglos = shAgglos.shAgglos;

% NOTE(amotta): This is necessary, so that we can mix the neurite
% agglomerates with the spine head agglomerates when building the LUTs.
shAgglos = cellfun(@double, shAgglos, 'UniformOutput', false);

%% Calculate positions
clear cur*;

curWmean = @(w, v) ...
    sum((w / sum(w, 1)) .* v, 1);
curShPos = round(cell2mat(cellfun( ...
    @(ids) curWmean(sizes(ids(:)), points(ids, :)), ...
    shAgglos, 'UniformOutput', false)));

curMask = isnan(adsT.pos);
curShPos = curShPos(adsT.shId, :);
adsT.pos(curMask) = curShPos(curMask);

%% Illustrate axon-spine interface
clear cur*;
curVxSize = param.raw.voxelSize;

for curIdx = 1:height(adsT)
    curAds = adsT(curIdx, :);
    
    %% Plot isosurface
    curAgglos = { ...
        conn.axons{curAds.axonId}; ...
        conn.dendrites{curAds.dendId}};
    
    curConfSegIds = sort(cell2mat(curAgglos(:)));
    curConfSegIds = curConfSegIds(diff(curConfSegIds) == 0);
    assert(isempty(curConfSegIds));
    
    curLUT = Agglo.buildLUT(maxSegId, curAgglos);
        
    curBox = 1E3 * curAds.isoFovUm ./ curVxSize;
    curBox = round(curAds.pos(:) + [-1, +1] / 2 .* curBox(:));
    
    curSeg = loadSegDataGlobal(param.seg, curBox);
    curSeg(curSeg ~= 0) = curLUT(curSeg(curSeg ~= 0));
    
    curIsos = cell(size(curAgglos));
    for curAggloIdx = 1:numel(curAgglos)
        curIso = curSeg == curAggloIdx;
        
        % See +connectEM/+Consistency/buildAxonSpineInterfaceAreas.m
        curIso = imclose(curIso, strel('cube', 3));
        curIso = isosurface(curIso, 0.5);
        
        curIso.vertices = curIso.vertices .* curVxSize;
        curIsos{curAggloIdx} = curIso;
    end
    
    curFig = figure();
    curAx = axes(curFig); %#ok
    hold(curAx, 'on');
    
    curPatches = cellfun( ...
        @(iso) patch(curAx, iso), curIsos);
    set(curPatches, ...
        'EdgeColor', 'none', ...
       {'FaceColor'}, num2cell(colors, 2));
    material(curPatches, 'dull');
    
    axis(curAx, 'equal');
    
    xlim(curAx, [0, diff(curBox(1, :)) * curVxSize(1)]);
    ylim(curAx, [0, diff(curBox(2, :)) * curVxSize(2)]);
    zlim(curAx, [0, diff(curBox(3, :)) * curVxSize(3)]);
    
    view(curAx, 3);
    camlight(curAx);
    
    curAx.XAxis.Visible = 'off';
    curAx.YAxis.Visible = 'off';
    curAx.ZAxis.Visible = 'off';
    
    curAx.BoxStyle = 'full';    
    curAx.Box = 'on';
    
    connectEM.Figure.config(curFig, info);
    
    %% Plot EM stack
    % See +connectEM/+Consistency/buildAxonSpineInterfaceAreas.m
    curAgglos = { ...
        conn.axons{curAds.axonId}; ...
        shAgglos{curAds.shId}};
    curLUT = Agglo.buildLUT(maxSegId, curAgglos);
        
    curBox = 1E3 * curAds.isoFovUm ./ curVxSize;
    curBox = round(curAds.pos(:) + [-1, +1] / 2 .* curBox(:));
    
    curSeg = loadSegDataGlobal(param.seg, curBox);
    curSeg(curSeg ~= 0) = curLUT(curSeg(curSeg ~= 0));
    
    curSeg(imclose(curSeg == 3, strel('cube', 3))) = 3;
    curSeg(imclose(curSeg == 2, strel('cube', 3))) = 2;
    curSeg(imclose(curSeg == 1, strel('cube', 3))) = 1;
    
    % NOTE(amotta): Closing doesn't work at end of cube.
    curSeg(1, :, :) = 3; curSeg(end, :, :) = 3;
    curSeg(:, 1, :) = 3; curSeg(:, end, :) = 3;
    curSeg(:, :, 1) = 3; curSeg(:, :, end) = 3;
    
   [curEdges, curBorder] = SynEM.Svg.findEdgesAndBorders(curSeg);
    curBorder = curBorder(all(curEdges == [1, 2], 2));
    assert(isscalar(curBorder));
    
    curBorderPos = curBorder.Centroid;
    curBorderVx = curBorder.PixelIdxList;
    
    curRaw = loadRawData(param.raw, curBox);
    curRawSize = size(curRaw);
    
    curRaw = double(curRaw);
    curRaw = max(curRaw - emRange(1), 0);
    curRaw = min(curRaw / diff(emRange), 1);
    
    curRaw = repmat(curRaw(:), [1, 3]);
    for curDimIdx = 1:3
        curRaw(curBorderVx, curDimIdx) = ...
            asiAlpha .* asiColor(curDimIdx) ...
          + (1 - asiAlpha) .* curRaw(curBorderVx, curDimIdx);
    end
    
    curRaw = reshape(curRaw, [curRawSize, 3]);
    curEmBox = 1E3 * curAds.emFovUm ./ curVxSize;
    curEmBox = round(curBorderPos(:) + [-1, +1] / 2 .* curEmBox(:));
    
    curRaw = curRaw( ...
        curEmBox(1, 1):curEmBox(1, 2), ...
        curEmBox(2, 1):curEmBox(2, 2), ...
        curEmBox(3, 1):curEmBox(3, 2), :);
    curRaw = permute(curRaw, [1, 2, 4, 3]);
    implay(curRaw);
end
