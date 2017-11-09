% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
axonFile = fullfile(rootDir, 'aggloState', 'axons_10_a.mat');
tempOutputDir = '/tmpscratch/amotta/l4/chiasma-detection';

chiasmaParam = struct;
chiasmaParam.minNrChiasmaExits = 3;
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

%% run chiasma detection
[job, runId] = connectEM.detectChiasmataSuperSuper( ...
    param, chiasmaParam, agglos, tempOutputDir);
fprintf('Running chiasma detection %s... ', runId);

wait(job);
fprintf('done!\n');

%% collecting results
fprintf('Collecting results... ')
chiasmata = connectEM.Chiasma.Detect.collectResults( ...
    tempOutputDir, runId, numel(agglos));
fprintf('done!\n');

%% evaluation
connectEM.Chiasma.Detect.evaluate(chiasmata, agglos);

%% building output
out = struct;
out.info = info;
out.runId = runId;
out.chiasmata = chiasmata;