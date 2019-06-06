% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
config = 'ex144_08x2_mrNet';

classEmParam = struct;
classEmParam.agglo = struct;
classEmParam.agglo.maxSegCount = 1;
classEmParam.agglo.padSize = [0, 0, 0];

fileName = 'segmentFeatures.mat';

info = Util.runInfo();
Util.showRunInfo(info);

config = loadConfig(config);
param = config.param;

%% Calculate features
job = TypeEM.Pipeline.buildFeatures(param, classEmParam, fileName);
Cluster.waitForJob(job);
