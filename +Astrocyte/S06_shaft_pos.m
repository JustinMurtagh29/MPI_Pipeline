rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

borderMeta = fullfile(rootDir, 'globalBorder.mat');
borderMeta = load(borderMeta, 'borderSize', 'borderCoM');
borderMeta = structfun(@double, borderMeta, 'UniformOutput', false);

load('/u/yyener/astrocyte/synapses/syn.mat');

edgeIds = syn.synapses.edgeIdx(syn.synapses.type == 'Shaft');

graphFile = fullfile(rootDir, 'graph.mat');
graph = load(graphFile, 'edges', 'borderIdx');


weightedMean = @(w, v) ...
    sum((w / sum(w, 1)) .* v, 1);

pos = cellfun( ...
    @(ids) weightedMean( ...
        borderMeta.borderSize(graph.borderIdx(ids)), ...
        borderMeta.borderCoM(graph.borderIdx(ids), :)), ...
	edgeIds, 'UniformOutput', false);
pos = round(cell2mat(pos));

%pos(cellfun(@isempty, borderIds), :) = nan;