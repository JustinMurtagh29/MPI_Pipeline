% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
outDir = '/tmpscratch/amotta/l4/2018-08-29-spine-head-isosurfaces';
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

% for speed
segParam = struct;
segParam.root = '/tmpscratch/amotta/l4/2012-09-28_ex145_07x2_ROI2017/segmentation/1';
segParam.backend = 'wkwrap';

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));

param = param.p;
param.seg = segParam;

maxSegId = Seg.Global.getMaxSegId(param);
conn = connectEM.Connectome.load(param, connFile);

shAgglos = load(shFile);
shAgglos = shAgglos.shAgglos;

heuristics = load(fullfile(rootDir, 'heuristicResult.mat'));
heuristics.segIds = double(heuristics.segIds);

%% Filtering
somaAgglos = conn.dendrites(conn.denMeta.targetClass == 'Somata');
nuclearSegIds = heuristics.segIds(heuristics.nucleiScore > 0.5);
vesselSegIds = heuristics.segIds(heuristics.vesselScore > 0.5);

keepMask = [somaAgglos; {nuclearSegIds; vesselSegIds}];
keepMask = not(Agglo.buildLUT(maxSegId, keepMask));

shAgglos = cellfun( ...
    @(segIds) segIds(keepMask(segIds)), ...
    shAgglos, 'UniformOutput', false);
shAgglos(cellfun(@isempty, shAgglos)) = [];

%% Build isosurfaces
Visualization.exportAggloToAmira( ...
    param, shAgglos, outDir, 'reduce', 0.05, ...
    'smoothSizeHalf', 4, 'smoothWidth', 8);
