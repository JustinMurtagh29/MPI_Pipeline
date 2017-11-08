% Collect results from chiasmata detection run `20171009T184902` and store
% it in a single file. This will speed up all further calculations.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
axonFile = fullfile(rootDir, 'aggloState', 'axons_06_c.mat');

chiasmaDir = '/tmpscratch/kboerg/chiasmata';
chiasmaId = '20171009T184902';

chiasmaParam = struct;
chiasmaParam.sphereRadiusOuter = 10000;
chiasmaParam.sphereRadiusInner = 1000;
chiasmaParam.minNodeDist = 2000;
chiasmaParam.clusterSize = 2000;

info = Util.runInfo();

%% load data
axonCount = load(axonFile, 'indBigAxons');
axonCount = sum(axonCount.indBigAxons);

%% collect chiasmata
chiasmata = ...
    connectEM.Chiasma.Detect.collectResults( ...
        chiasmaDir, chiasmaId, axonCount);

%% build output
out = struct;
out.info = info;
out.chiasmata = chiasmata;