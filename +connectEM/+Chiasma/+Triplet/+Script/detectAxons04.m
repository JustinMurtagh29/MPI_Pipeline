% For the published analysis of synapse size consistency, we've used an
% axon state that was split at all branch points and chiasmata. See:
%
%   connectEM.Chiasma.Triplet.Script.detectAxons19
%   git@gitlab.mpcdf.mpg.de:connectomics/pipeline.git 8aa31f9e06582be27b008d5d2dca3ceb70e3196b
%   amotta@gaba. MATLAB 9.3.0.713579 (R2017b). 04-Jun-2018 20:33:24
%
% Here, we try to replicate exactly this, but on the initial, fully-
% automatically generated axon agglomerates (i.e., "axons_04").
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
axonFile = fullfile(rootDir, 'aggloState', 'axons_04.mat');

outputDir = fullfile(rootDir, ...
    'tripletDetection', '20200608T164619-on-axons-04');
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

outFile = sprintf('%s_chiasmata.mat', runId);
outFile = fullfile(outputDir, outFile);

Util.saveStruct(outFile, out);
system(sprintf('chmod a-w "%s"', outFile));
