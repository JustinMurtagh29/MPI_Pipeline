% This script lists the isosurfaces of all axons that innervate somata.
% Because Heiko's isosurfaces are based on an older connectome, the indices
% need to be translated.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
plyDir = '/tmpscratch/amotta/l4/2018-01-24-axons-18a-isosurfaces/ply';
outFile = '/tmpscratch/amotta/l4/2018-01-24-axons-18a-isosurfaces/soma-targeting_v1.mat';

connOldFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a_ax_spine_syn_clust.mat');
connNewFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

connOld = load(connOldFile);
[connNew, synNew, axonClasses] = ...
    connectEM.Connectome.load(param, connNewFile);

maxSegId = Seg.Global.getMaxSegId(param);
segPoints = Seg.Global.getSegToPointMap(param);

%% Find soma-targeting axons
synNewT = connectEM.Connectome.buildSynapseTable(connNew, synNew);
synNewT.targetClass = connNew.denMeta.targetClass(synNewT.postAggloId);
synNewT(synNewT.targetClass ~= 'Somata', :) = [];

somaAxonIds = unique(synNewT.preAggloId);
inhAxonMask = connNew.axonMeta.axonClass(somaAxonIds) == 'Inhibitory';

%% Calculating synapse meta data
synNew.synapses.prePos = cellfun( ...
    @(segIds) mean(segPoints(segIds, :), 1), ...
    synNew.synapses.presynId, 'UniformOutput', false);
synNew.synapses.postPos = cellfun( ...
    @(segIds) mean(segPoints(segIds, :), 1), ...
    synNew.synapses.postsynId, 'UniformOutput', false);
synNew.synapses.pos = cell2mat(cellfun( ...
    @(prePos, postPos) (prePos + postPos) / 2, ...
    synNew.synapses.prePos, synNew.synapses.postPos, ...
    'UniformOutput', false));
synNewT.pos = synNew.synapses.pos(synNewT.id, :);
synNewT.isInh = ismember(synNewT.preAggloId, somaAxonIds(inhAxonMask));

%% Translate indices to old connectome
axonOldLUT = Agglo.buildLUT(maxSegId, connOld.axons);
somaAxonIds = cellfun( ...
    @(ids) mode(nonzeros(axonOldLUT(ids))), ...
    connOld.axons(somaAxonIds));

%% Build output
out = struct;
out.axonPlyFiles = arrayfun( ...
	@(id) sprintf('iso-%d.ply', id), ...
	somaAxonIds, 'UniformOutput', false);
out.axonPlyFiles = fullfile(plyDir, out.axonPlyFiles);
out.axonIsInh = inhAxonMask;

out.synapsePos = synNewT.pos;
out.synapseIsInh = synNewT.isInh;
out.info = info;

Util.saveStruct(outFile, out);
Util.protect(outFile);
