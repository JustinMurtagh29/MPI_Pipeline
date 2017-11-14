% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
axonFile = fullfile( ...
    rootDir, 'aggloState', ...
    'axons_11_a_plus_remaining_E3a.mat');

outputDir = fullfile( ...
    rootDir, 'chiasmataSplitting', ...
    '20171113T171102-on-axons-11a-plus-remaining-E3a');
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