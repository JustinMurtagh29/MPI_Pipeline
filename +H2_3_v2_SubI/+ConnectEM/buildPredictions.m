% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% HACKHACKHACK
% NOTE(amotta): This is a huge mess. The training data is located in my
% repository, SynEM is from Benedikt's repository, and the SynEM classifier
% is loaded from the pipeline repo.
%   Let's get rid of the dependency on the pipeline repo by building the
% feature map de-novo. The `paper` version is identical to the one stored
% in the pipeline repository up to different `border` values in the
% `AverageFilter`. But this doesn't affect `fm.invertDirection`.
addpath('/gaba/u/amotta/code/benedikt', '-end');
%{
%% Configuration
config = 'ex144_08x2_mrNet';
%}
rootDir = '/tmpscratch/sahilloo/data/H2_3_v2_U1_SubI/pipelineRun_mr2e_wsmrnet_HC/';
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

classifierFile = '/gaba/u/sahilloo/H2_3_v2_U1_SubI/connect-em/classifier.mat';

areaT = 10;
fmap = SynEM.getFeatureMap('paper');

info = Util.runInfo();
Util.showRunInfo(info);
%{
%% Loading data
config = loadConfig(config);
param = config.param;
%}
%% Build predictions
ConnectEM.Pipeline.buildPredictions( ...
    param, classifierFile, 'areaT', areaT, 'fmap', fmap);
