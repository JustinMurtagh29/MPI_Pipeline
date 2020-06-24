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
classEmParam.agglo.maxSegCount = 1;
classEmParam.agglo.padSize = [0, 0, 0];

% introduced later by AM for mSEM - try for SBEM
classEmParam.texture.inputs = struct;
classEmParam.texture.inputs.name = 'raw'; % no 'membrane' based affinity features for SBEM
classEmParam.texture.inputs.load = @(p, box) ...
        (single(wkwLoadRoi(p.raw.root, box)) - p.norm.meanVal) / p.norm.stdVal;

fileName = 'segmentFeatures.mat';

info = Util.runInfo();
Util.showRunInfo(info);

%% Calculate features
job = TypeEM.Pipeline.buildFeatures(param, classEmParam, fileName);
Cluster.waitForJob(job);
