% Run chiasma detection on a random selection of agglomerates. This script
% is intended for debugging and / or benchmarking purposes.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
axonFile = fullfile(rootDir, 'aggloState', 'axons_10_a_plus_10kE3a_b.mat');

outputDir = fullfile( ...
    rootDir, 'chiasmataSplitting', ...
    '20171104T181213-on-axons-10a-plus-10kE3a');
tempOutputDir = '/tmpscratch/amotta/l4/chiasma-detection';

chiasmaParam = struct;
chiasmaParam.sphereRadiusOuter = 10000;
chiasmaParam.sphereRadiusInner = 1000;
chiasmaParam.minNodeDist = 2000;
chiasmaParam.clusterSize = 2000;

info = Util.runInfo();

%% load data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

agglos = load(axonFile, 'axons', 'indBigAxons');
agglos = agglos.axons(agglos.indBigAxons);

%% prepare parameter
chiasmaParam.voxelSize = param.raw.voxelSize;

chiasmaParamKvPairs = cat( ...
    2, fieldnames(chiasmaParam), struct2cell(chiasmaParam));
chiasmaParamKvPairs = transpose(chiasmaParamKvPairs);

p = param;
p = Util.modifyStruct(p, chiasmaParamKvPairs{:});

%% run on
nodeCount = arrayfun(@(a) size(a.nodes, 1), agglos);
[~, aggloIds] = sort(nodeCount, 'descend');
aggloIds = aggloIds(1);

tic;
profile on;
outputs = arrayfun(@(idx) connectEM.detectChiasmata( ...
    p, agglos(idx).nodes(:, 1:3), agglos(idx).edges, false, []), aggloIds);
toc;