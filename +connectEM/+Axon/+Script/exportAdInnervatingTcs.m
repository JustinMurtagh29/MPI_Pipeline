% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');
outputDir = '/home/amotta/Desktop/ad-innervating-tc-axons';

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
clear cur*;
synT = connectEM.Connectome.buildSynapseTable(conn, syn);
synT.axonClass = conn.axonMeta.axonClass(synT.preAggloId);
synT.targetClass = conn.denMeta.targetClass(synT.postAggloId);

synT = synT( ...
    synT.axonClass == 'Thalamocortical' ...
  & synT.targetClass == 'ApicalDendrite', :);

tcAdT = table;
[tcAdT.axonId, ~, curUniIds] = unique(synT.preAggloId);
tcAdT.synIds = accumarray(curUniIds, synT.id, [], @(ids) {ids});
tcAdT.dendIds = accumarray(curUniIds, synT.postAggloId, [], @(ids) {unique(ids)});

%% Generate NML files
clear cur*;
mkdir(outputDir);

% Export random examples
rng(0);
curExportIds = randperm(height(tcAdT));
curExportIds = curExportIds(1:50);

curDigits = ceil(log10(1 + numel(curExportIds)));

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

for curIdx = 1:numel(curExportIds)
    curRow = curExportIds(curIdx);
    curAxonId = tcAdT.axonId(curRow);
    curDendIds = tcAdT.dendIds{curRow};    
    curSynIds = tcAdT.synIds{curRow};
    
    curAgglos = cellfun(@vertcat, ...
        syn.synapses.presynId(curSynIds), ...
        syn.synapses.postsynId(curSynIds), ...
        'UniformOutput', false);
    curAgglos = [ ...
        conn.axons(curAxonId); ...
        conn.dendrites(curDendIds); ...
        curAgglos(:)];
    
    curSkel = skel;
    curSkel = Skeleton.fromMST(cellfun( ...
        @(ids) segPoints(ids, :), curAgglos, ...
        'UniformOutput', false), param.raw.voxelSize, curSkel);
    
    curSkel.names{1} = sprintf('Axon %d', curAxonId);
    curSkel.names(1 + (1:numel(curDendIds))) = arrayfun( ...
        @(id) sprintf('%d (%s)', id, conn.denMeta.targetClass(id)), ...
        curDendIds, 'UniformOutput', false);
        
    curSkel.names((2 + numel(curDendIds)):end) = arrayfun( ...
        @(id) sprintf('Synapse %d', id), curSynIds, 'UniformOutput', false);
    
    curSkelFile = sprintf( ...
        '%0*d_axon-%d.nml', curDigits, curIdx, curAxonId);
    curSkel.write(fullfile(outputDir, curSkelFile));
end
