% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/tmpscratch/amotta/l23/2018-10-09-mrnet-pipeline-run';

param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

margin = param.tileBorder(:, 2);
assert(isequal(margin, [256; 256; 128]));

info = Util.runInfo();
Util.showRunInfo(info);

%% Run per-cube agglomeration
taskArgs = arrayfun( ...
    @(local) {local.bboxSmall}, ...
    param.local, 'UniformOutput', false);
sharedArgs = {param, 'margin', margin};
sharedArgLocs = [1, 3, 4];

job = Cluster.startJob( ...
    @connectEM.Hierarchical.runBox, ...
    taskArgs, 'numOutputs', 2, ...
    'sharedInputs', sharedArgs, ...
    'sharedInputsLocation', sharedArgLocs, ...
    'cluster', {'priority', 100, 'time', '1:00:00', 'memory', 48});
Cluster.waitForJob(job);
