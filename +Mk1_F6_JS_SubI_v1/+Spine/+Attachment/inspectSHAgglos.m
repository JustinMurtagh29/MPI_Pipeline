% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
% Modified by
%   Sahil Loomba <sahil.loomba@brain.mpg.de>

clear;

%% configuration
rootDir = '/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet/';
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

nmlDir = fullfile(param.saveFolder, 'tracings', 'spine-attachment');
dendFile = fullfile(param.saveFolder, 'aggloMat', '20191227T134319_ga_20191224T235355optimParams_agglomeration/20191227T220548_results.mat');

parameters.scale.x = num2str(param.raw.voxelSize(1));
parameters.scale.y = num2str(param.raw.voxelSize(2));
parameters.scale.z = num2str(param.raw.voxelSize(3));
parameters.offset.x = '0';
parameters.offset.y = '0';
parameters.offset.z = '0';
parameters.experiment.name = param.experimentName;

aggloParam = struct;
% NOTE(amotta): The spine head probability was optimized using:
% L23.TypeEM.buildSegmentClassifier
% git@gitlab.mpcdf.mpg.de:connectomics/amotta.git 6a46f4500a4392aab1a5c33b8f0b4f496d1f22b5 (dirty)
% amotta@gaba01. MATLAB 9.3.0.713579 (R2017b). 11-Dec-2019 15:19:38
%
% More detailed description in
% https://mhlablog.net/2020/06/09/nhp-spine-head-detection/
aggloParam.minSpineHeadProb = 0.661619; % 0.5545;
aggloParam.maxVesselScore = 0.5;
aggloParam.maxNucleusScore = 0.5;

aggloParam.minEdgeProb = 0.98;

outputFolder = fullfile(param.saveFolder, 'spineAttachment', ...
        sprintf('spineheadProb-%.02f-edgeProb-%.02f', aggloParam.minSpineHeadProb, aggloParam.minEdgeProb));

info = Util.runInfo();

Util.log('Reading dendrites from "%s"', dendFile);

%% actually do some work
Util.log('Building spine heads');
sh = struct;
[sh.agglos, sh.edges] = ...
    L4.Spine.Head.buildAgglos(param, aggloParam);

Util.log('Loading graph');
graph = fullfile(rootDir, 'graph.mat');
graph = load(graph, 'edges', 'prob');
segmentMeta = load([param.saveFolder 'segmentMeta.mat'], 'voxelCount', 'point', 'maxSegId');
segmentMeta.point = segmentMeta.point';

Util.log('Export skeletons:')
agglosOut = sh.agglos(1:100);
Superagglos.skeletonFromAgglo(graph.edges, segmentMeta, ...
    agglosOut, 'shAgglos', outputFolder, parameters);

Util.log('Done');
