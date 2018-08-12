% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');
outputDir = '/home/amotta/Desktop/AIS-innervating-axons';

minSynPre = 10;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

segPoints = Seg.Global.getSegToPointMap(param);

[conn, syn, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);

[conn, axonClasses] = ...
    connectEM.Connectome.prepareForSpecificityAnalysis( ...
        conn, axonClasses, 'minSynPre', minSynPre);

%% Find axons with AIS synapses
synT = connectEM.Connectome.buildSynapseTable(conn, syn);
synT.targetClass = conn.denMeta.targetClass(synT.postAggloId);
synT = synT(synT.targetClass == 'AxonInitialSegment', :);

aisT = table;
[aisT.axonId, ~, synIds] = unique(synT.preAggloId);
aisT.synIds = accumarray(synIds, synT.id, [], @(ids) {ids});

%% Generate NML files
clear cur*;
mkdir(outputDir);

% Export random examples
rng(0);
curExportIds = randperm(height(aisT));
curExportIds = curExportIds(1:50);

%{
% Export by decreasing number of AIS synapses
[curAisSynCounts, curExportIds] = sort( ...
    cellfun(@numel, aisT.synIds), 'descend');
curExportIds = curExportIds(curAisSynCounts > 1);
%}

curDigits = ceil(log10(1 + numel(curExportIds)));

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

for curIdx = 1:5%numel(curExportIds)
    curAisRow = curExportIds(curIdx);
    curAxonId = aisT.axonId(curAisRow);
    curSynIds = aisT.synIds{curAisRow};
    
    curAgglos = cellfun(@vertcat, ...
        syn.synapses.presynId(curSynIds), ...
        syn.synapses.postsynId(curSynIds), ...
        'UniformOutput', false);
    curAgglos = [conn.axons(curAxonId); curAgglos(:)];
    
    curSkel = skel;
    curSkel = Skeleton.fromMST(cellfun( ...
        @(ids) segPoints(ids, :), curAgglos, ...
        'UniformOutput', false), param.raw.voxelSize, curSkel);
    
    curSkel.names{1} = sprintf('Axon %d', curAxonId);
    curSkel.names(2:end) = arrayfun(@(id) ...
        sprintf('Synapse %d', id), curSynIds, ...
        'UniformOutput', false);
    
    curSkelFile = sprintf( ...
        '%0*d_ais-axon-%d.nml', curDigits, curIdx, curAxonId);
    curSkel.write(fullfile(outputDir, curSkelFile));
end
