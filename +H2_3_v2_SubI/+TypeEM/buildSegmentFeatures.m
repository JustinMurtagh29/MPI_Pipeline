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
classEmParam.agglo.maxSegCount = 1;
classEmParam.agglo.padSize = [0, 0, 0];

fileName = 'segmentFeatures.mat';

info = Util.runInfo();
Util.showRunInfo(info);

%% Calculate features
job = TypeEM.Pipeline.buildFeatures(param, classEmParam, fileName);
Cluster.waitForJob(job);
