% Calculate features for cubes that are used in the GT while training

rootDir = ['/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet/'];
load(fullfile(rootDir,'allParameter.mat'))
param = p;
param.experimentName = 'Mk1_F6_JS_SubI_v1_mrnet_wsmrnet';

classEmParam = struct;
classEmParam.agglo = struct;
classEmParam.agglo.padSize = [256, 256, 128];
classEmParam.agglo.minEdgeProb = 0.8; % 0.5
classEmParam.agglo.maxSegCount = 5;

fileName = 'segmentAgglomerateV2Features.mat'; %'segmentAgglomerateFeatures.mat';

% load cubeIds.mat from disk
nmlDir = fullfile(param.saveFolder, ...
     'tracings', 'box-seeded','spine-head-ground-truth', 'v4');
m = load(fullfile(nmlDir, 'cubeIds.mat'));
cubeIds = m.cubeIds;

% build job input arguments for GT cubeIds
taskSeparateInputs = arrayfun(@(i) {{i}}, cubeIds);
taskSharedInputs = {param, classEmParam, fileName};

% build and launch job
% * for single segments:
%   successfully used < 18 GB + ~3 hours
%   good evidence that < 12 GB + ~3 hours could work
% * for five segment agglomerates
%   ~42 GB + ~4:15 h
job = Cluster.startJob( ...
    @main, taskSeparateInputs, ...
    'sharedInputs', taskSharedInputs, ...
    'cluster', {'memory', 48, 'time', '24:00:00', 'priority',100}, ...
    'name', mfilename());

function main(param, classEmParam, fileName, cubeIdx)
    info = Util.runInfo();
    info.param = rmfield(info.param, 'param');
    
    param.classem = classEmParam;
    cubeParam = param.local(cubeIdx);
    rootDir = cubeParam.saveFolder;
    box = cubeParam.bboxSmall;
    
    % calculate features
    out = TypeEM.buildFeatures(param, box);
    out.gitHash = Util.getGitHash();
    out.info = info;
    
    % store result
    if ~exist(rootDir, 'dir'); mkdir(rootDir); end
    outFile = fullfile(rootDir, fileName);
    Util.saveStruct(outFile, out);
    Util.protect(outFile);
end

