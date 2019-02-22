% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
outputMapFile = '/tmpscratch/amotta/l4/2018-07-26-tracing-based-output-maps/20190117T143833_results.mat';

l4ConnRunId = '20190221T112510';
outDir = '/home/amotta/Desktop';

checkedSynIds = [ ...
    ... All primary spine synapses from L4 → AD (19.02.2019)
    7975, 68992, 77903, 185118, 186018, 97233, 119754, ...
    122684, 271912, 300946, 344896, 294170, 140701, 202773, ...
    ... Randomly selected primary spine synapses from L4 → L4 (20.02.2019)
    97662, 109112, 72299, 101803, 228602, 17999, ...
    275717, 163160, 44870, 163175, 33945, 94428];

info = Util.runInfo();
Util.showRunInfo(info);

%% Load data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

segPoints = Seg.Global.getSegToPointMap(param);
segSizes = Seg.Global.getSegToSizeMap(param);

curData = load(outputMapFile);
axonData = curData.axonData;

conn = curData.info.param.connFile;
conn = connectEM.Connectome.load(param, conn);

[curDir, curFile] = fileparts(outputMapFile);
curAsiFile = sprintf('%s__%s_connectome.mat', curFile, l4ConnRunId);
curData = load(fullfile(curDir, curAsiFile));

l4SynT = curData.synT;
l4AsiT = curData.asiT;

l4AsiT = l4AsiT(l4AsiT.area > 0, :);
l4AsiT = connectEM.Consistency.Calibration.apply(l4AsiT);

%% Prepare synapse table
l4SynT.targetClass = conn.denMeta.targetClass(l4SynT.postAggloId);
[~, curIds] = ismember(l4AsiT.id, l4SynT.id);

l4SynT.type(:) = categorical({'Shaft'});
l4SynT.type(curIds) = l4AsiT.type;

l4SynT.asiArea(:) = nan;
l4SynT.asiArea(curIds) = l4AsiT.area;

%% Select random synapses to query
clear cur*;
rng(0);

curSynT = l4SynT;
curSynT = curSynT(curSynT.type == 'PrimarySpine', :);

% 20.02.2019: All L4 → AD and 20 L4 → WC synapses
querySynIds = checkedSynIds;
curSynT(ismember(curSynT.id, querySynIds), :) = [];

% 22.02.2019: Select 40 more L4 → WC synapses
curSynIds = curSynT.id(curSynT.targetClass == 'WholeCell');
curSynIds = curSynIds(randperm(numel(curSynIds)));
curSynIds = curSynIds(1:40);

querySynIds = union(querySynIds, curSynIds, 'stable');
curSynT(ismember(curSynT.id, querySynIds), :) = [];

%% Generate table for Heiko
clear cur*;

queryT = table;
queryT.synId = querySynIds(:);

[~, curIds] = ismember(queryT.synId, l4SynT.id);
[~, queryT.nmlFile] = arrayfun( ...
    @(axonId) fileparts(axonData(axonId).nmlFile), ...
    l4SynT.preAggloId(curIds), 'UniformOutput', false);

format long;
fprintf('Queries\n\n');
disp(queryT);
