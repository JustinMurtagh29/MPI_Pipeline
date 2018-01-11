% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
cacheDir = '/tmpscratch/amotta/l4/2018-01-axon-16b-flight-paths';
axonFile = fullfile(rootDir, 'aggloState', 'axons_16_b.mat');

% For now, let's use the same neighborhood size as in
% `+connectEM/axonAggloStateVisualization.m`.
nhoodSize = 3;

[~, axonName] = fileparts(axonFile);
cacheFile = fullfile(cacheDir, strcat(axonName, '_flights.mat'));
clear axonName;

info = Util.runInfo();

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

axons = load(axonFile);

%% use WKW segmentation
param.seg = struct;
param.seg.root = '/tmpscratch/amotta/l4/2012-09-28_ex145_07x2_ROI2017/segmentation/1';
param.seg.backend = 'wkwrap';

%% look up flight paths
[axonFlights, axonFlightsMeta] = ...
    Superagglos.getFlightPathSegIds( ...
        param, axons.axons(axons.indBigAxons), nhoodSize);
save(cacheFile, '-v7.3', '-append', 'axonFlights', 'axonFlightsMeta');