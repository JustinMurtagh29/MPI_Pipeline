% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
outDir = '/tmpscratch/amotta/l4/2018-08-29-spine-head-isosurfaces';
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

numGroups = 100;

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

%% Group spine heads for performance reasons
rng(0);
randIds = randperm(numel(shAgglos));

shGroups = ceil(linspace( ...
    1, numel(shAgglos), numGroups + 1));
shGroups = arrayfun( ...
    @(a, b) cell2mat(shAgglos(randIds(a:b))), ...
    shGroups(1:(end - 1)), shGroups(2:end), ...
    'UniformOutput', false);

%% Build isosurfaces
Visualization.exportAggloToAmira( ...
    param, shGroups, outDir, 'reduce', 0.05, ...
    'smoothSizeHalf', 4, 'smoothWidth', 8);
