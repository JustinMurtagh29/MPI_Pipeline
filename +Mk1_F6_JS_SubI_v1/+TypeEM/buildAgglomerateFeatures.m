% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = ['/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet/'];
load(fullfile(rootDir,'allParameter.mat'))
param = p;
param.experimentName = 'Mk1_F6_JS_SubI_v1_mrnet_wsmrnet';

classEmParam = struct;
classEmParam.agglo = struct;
classEmParam.agglo.padSize = [256, 256, 128];
classEmParam.agglo.minEdgeProb = 0.5;
classEmParam.agglo.maxSegCount = 5;

fileName = 'segmentAgglomerateFeatures.mat';

info = Util.runInfo();
Util.showRunInfo(info);

config = loadConfig(config);
param = config.param;

%% Calculate features
job = TypeEM.Pipeline.buildFeatures(param, classEmParam, fileName);
Cluster.waitForJob(job);
