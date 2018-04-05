% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
synFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v3_ax_spine_clustered.mat');
connName = 'connectome_axons_18_a_ax_spine_syn_clust';
outDir = '/home/amotta/Desktop/soma-specific-axons';

minSynPre = 10;
minSomaSpec = 0.3;

% Set to export only these axons
chosenAxonIds = [630]; %#ok

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

syn = load(synFile);
conn = connectEM.Connectome.load(param, connName);
segPoints = Seg.Global.getSegToPointMap(param);

%% Identify inhibitory axons
axonClasses = ...
    connectEM.Connectome.buildAxonClasses(conn, 'minSynPre', minSynPre);
inhAxonIds = axonClasses(2).axonIds;

%% Identify soma-specific axons
[classConnectome, targetClasses] = ...
    connectEM.Connectome.buildClassConnectome(conn);
somaMask = (targetClasses == 'Somata');

inhSomaSpecs = classConnectome(inhAxonIds, :);
inhSomaSpecs = inhSomaSpecs(:, somaMask) ./ sum(inhSomaSpecs, 2);

inhSomaSpecAxonIds = inhAxonIds;
inhSomaSpecAxonIds(inhSomaSpecs < minSomaSpec) = [];

%% Export to webKNOSSOS
if ~exist('chosenAxonIds', 'var') ...
        || isempty(chosenAxonIds)
	exportAxonIds = inhSomaSpecAxonIds;
else
    exportAxonIds = chosenAxonIds(:);
end

synT = connectEM.Connectome.buildSynapseTable(conn, syn);
synT.ontoSoma = conn.denMeta.targetClass(synT.postAggloId);
synT.ontoSoma = (synT.ontoSoma == 'Somata');

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

somaStrings = {'', ' (Onto soma)'};
numDigits = ceil(log10(1 + numel(exportAxonIds)));

for curIdx = 1:numel(exportAxonIds)
    curId = exportAxonIds(curIdx);
    curAxon = conn.axons(curId);
    
    curSynT = synT(synT.preAggloId == curId, :);
    curSomaIds = unique(curSynT.postAggloId(curSynT.ontoSoma));
    
    curSyn = cellfun(@vertcat, ...
        syn.synapses.presynId(curSynT.id), ...
        syn.synapses.postsynId(curSynT.id), ...
        'UniformOutput', false);
    curSomata = conn.dendrites(curSomaIds);
    
    curNodes = cellfun( ...
        @(segIds) segPoints(segIds, :), ...
        vertcat(curAxon, curSomata), ...
        'UniformOutput', false);
    curNames = [ ...
       {sprintf('Axon %d', curId)}; ...
        arrayfun( ...
            @(id) sprintf('Soma (Dendrite %d)', id), ...
            curSomaIds, 'UniformOutput', false); ...
        arrayfun( ...
            @(id, somaFlag) sprintf( ...
                'Synapse %d%s', id, somaStrings{1 + somaFlag}), ...
            curSynT.id, curSynT.ontoSoma, 'UniformOutput', false)];
    
    curSkel = Skeleton.fromMST( ...
        curNodes, param.raw.voxelSize, skel);
    curSkel = curSkel.addNodesAsTrees( ...
        cell2mat(cellfun(@(segIds) ...
            mean(segPoints(segIds, :), 1), ...
            curSyn, 'UniformOutput', false)));
    curSkel.names = curNames;
    
    curSkelName = sprintf('%0*d_axon-%d.nml', numDigits, curIdx, curId);
    curSkel.write(fullfile(outDir, curSkelName));
end
