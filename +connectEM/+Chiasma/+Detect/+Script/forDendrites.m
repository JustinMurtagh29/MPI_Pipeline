% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
dendriteFile = fullfile(rootDir, 'aggloState', 'dendrites_flight_02.mat');

tempOutputDir = '/tmpscratch/scchr/L4/chiasmata/dendrites/round5/';
outputDir = '/tmpscratch/scchr/L4/chiasmata/dendrites/evaluation/round5/';

chiasmaParam = struct;
chiasmaParam.sphereRadiusOuter = 22000;
chiasmaParam.sphereRadiusInner = 2000;
chiasmaParam.minNodeDist = 4000;
chiasmaParam.clusterSize = 2000;

info = Util.runInfo();

%% load data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

agglos = load(dendriteFile, 'dendrites', 'indBigDends');
agglos = agglos.dendrites(agglos.indBigDends);

% Edges = structfun(@(x)AxonEndings.minimalSpanningTree(x.nodes(:,1:3)),agglos);
Edges = [];
for i=1:length(agglos)
    agglos(i).edges = AxonEndings.minimalSpanningTree(agglos(i).nodes(:,1:3));
end
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
