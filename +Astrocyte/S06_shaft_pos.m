rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

borderMeta = fullfile(rootDir, 'globalBorder.mat');
borderMeta = load(borderMeta, 'borderSize', 'borderCoM');
borderMeta = structfun(@double, borderMeta, 'UniformOutput', false);

borderIds = syn.synapses.edgeIdx(syn.synapses.type == 'Shaft');

weightedMean = @(w, v) ...
    sum((w / sum(w, 1)) .* v, 1);

pos = cellfun( ...
    @(ids) weightedMean( ...
        borderMeta.borderSize(ids), ...
        borderMeta.borderCoM(ids, :)), ...
	borderIds, 'UniformOutput', false);
pos = round(cell2mat(pos));

pos(cellfun(@isempty, borderIds), :) = nan;