% NOTE: This script was written in such a way that additional query rounds
% can easily be implemented by replacing the paths in axonNmlFiles with the
% results from the previous annotation round. Slight changes to the query
% selection section should be all that's neeeded to get this going.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
% NOTE(amotta): This file is identical to 20180726T190355_results.mat with
% the exception of the additional `outputMap.axonData.segIds` field.
outputConnFile = '/tmpscratch/amotta/l4/2018-07-26-tracing-based-output-maps/20190117T143833_results.mat';
outputConnRunId = '20190221T112510';

queryCount = 50;
queryDir = '/home/amotta/Desktop';

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

segPoints = Seg.Global.getSegToPointMap(param);

% Loading tracings and synapses of L4 cells
[curDir, curFile] = fileparts(outputConnFile);
curFile = sprintf('%s__%s_connectome.mat', curFile, outputConnRunId);
curFile = fullfile(curDir, curFile);

curData = load(curFile);
synT = curData.synT;

outputMapFile = curData.info.param.outputMapFile;
curData = load(outputMapFile);
axonData = curData.axonData;
axonNmlFiles = {axonData.nmlFile};

connFile = curData.info.param.connFile;
[~, synAgglos] = connectEM.Connectome.load(param, connFile);

synAgglos = cellfun(@union, ...
    synAgglos.synapses.presynId, ...
    synAgglos.synapses.postsynId, ...
    'UniformOutput', false);

%% Select synapses to query
clear cur*;
rng(0);

queries = struct;
queries.spine = find(synT.isSpine);
queries.shaft = find(not(synT.isSpine));

% Shuffle
queries = structfun( ...
    @(v) v(randperm(numel(v))), ...
    queries, 'UniformOutput', false);

% Select subset
queries = structfun( ...
    @(v) v(1:min(numel(v), queryCount)), ...
    queries, 'UniformOutput', false);

%% Generate NML files
clear cur*;
curSynClasses = {'Shaft', 'Spine'};

for curAxonId = 1:numel(axonNmlFiles)
    curAxonNmlFile = axonNmlFiles{curAxonId};
   [~, curAxonName] = fileparts(curAxonNmlFile);
    
    curSynT = queries;
    curSynT = structfun( ...
        @(v) v(synT.preAggloId(v) == curAxonId), ...
        curSynT, 'UniformOutput', false);
    
    curSynT = cell2mat(struct2cell(curSynT));
    curSynT = synT(curSynT, :);
    
    curSkel = skeleton(curAxonNmlFile);
    curSkel = curSkel.deleteTreeWithName('Dendrite');
    
    for curSynIdx = 1:height(curSynT)
        curSyn = curSynT(curSynIdx, :);
        curSegIds = synAgglos{curSyn.id};
        
        curSkel = Skeleton.fromMST( ...
            segPoints(curSegIds, :), ...
            param.raw.voxelSize, curSkel);
        curSkel.names{end} = sprintf( ...
            'Synapse %d. %s. Todo', curSyn.id, ...
            curSynClasses{1 + curSyn.isSpine});
    end
    
    curSkel = Skeleton.setDescriptionFromRunInfo(curSkel, info);
    curSkel.write(fullfile(queryDir, sprintf('%s.nml', curAxonName)));
end
