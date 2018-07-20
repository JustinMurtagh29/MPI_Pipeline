% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');
synFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v3_ax_spine_clustered_classified.mat');
outDir = '/tmpscratch/amotta/l4/2018-07-17-axons-19a-without-soma-overlaps';

isoDir = [];
axonClassFile = fullfile(outDir, 'axonClasses_v4.mat');

minBorderDistUm = 5;

info = Util.runInfo();

%% load parameters
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

box = param.bbox;
voxelSize = param.raw.voxelSize;
maxSegId = Seg.Global.getMaxSegId(param);
segPoints = Seg.Global.getSegToPointMap(param);

[conn, ~, axonClasses] = ...
    connectEM.Connectome.load(param, connFile, synFile);
axonClasses = rmfield(axonClasses, 'nullAxonIds');

% for speed
param.seg = struct;
param.seg.root = '/tmpscratch/amotta/l4/2012-09-28_ex145_07x2_ROI2017/segmentation/1';
param.seg.backend = 'wkwrap';

%% prepare axon agglomerates
axons = conn.axons;

% Remove soma segments from axons
somaMask = conn.dendrites(conn.denMeta.targetClass == 'Somata');
somaMask = logical(Agglo.buildLUT(maxSegId, somaMask));

axons = cellfun( ...
    @(segIds) segIds(~somaMask(segIds)), ...
    axons, 'UniformOutput', false);

%% build class of axons that reach a minimum distance from dataset border
boxCore = round(1E3 * minBorderDistUm ./ param.raw.voxelSize);
boxCore = transpose(box) - [-1; +1] .* boxCore;

pointsInCore = @(pts) all(pts > boxCore(1, :) & pts < boxCore(2, :), 2);
aggloInCore = @(segIds) any(pointsInCore(segPoints(segIds, :)));

coreAxonClass = struct;
coreAxonClass.axonIds = find(cellfun(aggloInCore, axons));
coreAxonClass.title = 'core axons';

axonClasses(end + 1) = coreAxonClass;

%% generate isosurfaces
if ~isempty(isoDir)
    Visualization.exportAggloToAmira( ...
        param, axons, isoDir, 'reduce', 0.05, ...
        'smoothSizeHalf', 4, 'smoothWidth', 8);
end

%% generate axon classes
if ~isempty(axonClassFile)
    % Build paths
    plyFiles = arrayfun( ...
        @(axonId) sprintf('iso-%d.ply', axonId), ...
        transpose(1:numel(axons)), 'UniformOutput', false);
    plyFiles = fullfile(fileparts(axonClassFile), 'ply', plyFiles);
    
    % Build output
    out = struct;
    out.info = info;
    
    for axonClass = reshape(axonClasses, 1, [])
        axonClassName = strsplit(axonClass.title);
        axonClassName = axonClassName{1};
        
        out.(axonClassName) = plyFiles(axonClass.axonIds);
    end
    
    Util.saveStruct(axonClassFile, out);
    Util.protect(axonClassFile);
end
