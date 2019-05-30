% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');
outDir = '/tmpscratch/amotta/l4/2019-05-30-synapse-table-for-laura';

runId = datestr(now, 30);
info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, syn, axonClasses] = connectEM.Connectome.load(param, connFile);

%% Building synapse table
synT = connectEM.Connectome.buildSynapseTable(conn, syn);

% Synapse type
synT(:, 'isSpine') = [];
synT.type = syn.synapses.type(synT.id);

% Axon class
synT.axonClass = conn.axonMeta.axonClass(synT.preAggloId);

% Position
calculatePositions = @(n) feval(@(p) p(synT.id, :), ...
    connectEM.Synapse.calculatePositions(param, syn, n));
synT.posPre = calculatePositions('pre');
synT.posSyn = calculatePositions('prePost');
synT.posPost = calculatePositions('post');

%% Save output
out = struct;
out.info = info;
out.synapseTable = synT;

outFile = sprintf('%s_synapse-table.mat', runId);
outFile = fullfile(outDir, outFile);

Util.saveStruct(outFile, out);
Util.protect(outFile);
