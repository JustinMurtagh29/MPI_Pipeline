% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');
outDir = '/tmpscratch/amotta/l4/2019-11-06-connectome-for-specificity-movie';

minSynPre = 10;

info = Util.runInfo();
Util.showRunInfo(info);

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, syn, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);
[conn, axonClasses] = ...
    connectEM.Connectome.prepareForSpecificityAnalysis( ...
        conn, axonClasses, 'minSynPre', minSynPre);

%% Build synapse table
% As in /tmpscratch/amotta/l4/2019-07-29-same-axon-same-dendrite-pair-table/20190729T143522_sasd-pair-table.mat
synT = connectEM.Connectome.buildSynapseTable(conn, syn);
synT.isSpine = [];

synT = cat(2, synT, syn.synapses(synT.id, 'type'));
synT = synT(:, [1, end, 2:(end - 1)]);

synT.shId(:) = nan; % TODO(amotta)
synT.synIds = num2cell(synT.id);
synT.borderIds = syn.synapses.edgeIdx(synT.id);

synT.axonClass = conn.axonMeta.axonClass(synT.preAggloId);
synT.targetClass = conn.denMeta.targetClass(synT.postAggloId);

curPos = connectEM.Synapse.calculatePositions(param, syn);
synT.pos = round(curPos(synT.id, :));

synT.area(:) = nan; % TODO(amotta)
synT.postCellId = conn.denMeta.cellId(synT.postAggloId);

%% Save output
clear cur*

outFile = fullfile(outDir, sprintf( ...
    '%s_synapse-table.mat', datestr(now, 30)));

out = struct;
out.info = info;
out.synT = synT;

Util.saveStruct(outFile, out);
Util.protect(outFile);
