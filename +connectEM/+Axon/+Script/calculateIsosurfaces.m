% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');
synFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v3_ax_spine_clustered_classified.mat');

isoDir = '/tmpscratch/amotta/l4/2018-07-17-axons-19a-without-soma-overlaps';
axonClassFile = '';

info = Util.runInfo();

%% load parameters
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);

[conn, ~, axonClasses] = ...
    connectEM.Connectome.load(param, connFile, synFile);

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
