% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/tmpscratch/sahilloo/data/H2_3_v2_U1_SubI/pipelineRun_mr2e_wsmrnet/';
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;
param.experimentName = 'H2_3_v2_U1_SubI_mr2e_wsmrnet';

classEmParam = struct;
classEmParam.agglo = struct;
classEmParam.agglo.padSize = [256, 256, 128];
classEmParam.agglo.minEdgeProb = 0.5;
classEmParam.agglo.maxSegCount = 5;

fileName = 'segmentAgglomerateFeatures.mat';

info = Util.runInfo();
Util.showRunInfo(info);

param = config.param;

%% Calculate features
job = TypeEM.Pipeline.buildFeatures(param, classEmParam, fileName);
Cluster.waitForJob(job);
