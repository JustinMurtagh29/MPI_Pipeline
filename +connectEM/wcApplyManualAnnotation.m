function [ wholeCells,dendrites,indBigDends,myelinDend ] = wcApplyManualAnnotation(dendritesAndWC, correctionFolder, axons,somaSegIds,somaLUT,sMpoints )
% This function does the rounds of manual corrections (splitting and
% extending) for the whole cells

% the corrections are based on checking the output of
% connectEM.searchSkippedAgglos


[dendrites,dendriteLUT] = connectEM.applyAggloCorrections(dendritesAndWC,p,correctionFolder,0,sMpoints,axons,1);

% test for duplets
agglos = cell2mat(Superagglos.transformAggloNewOldRepr(dendrites));
assert(numel(agglos)==numel(unique(agglos)))

dendriteSegIds = find(dendriteLUT);
[ismem,ind] = ismember(somaSegIds,dendriteSegIds);
% get each dend id which contains most of the seg ids of each soma, in
% this last step do not make unique, as all cells are now supposed to
% not be merged and thus no duplets exist
wholeCellId = (accumarray(somaLUT(somaSegIds(ismem))',dendriteLUT(dendriteSegIds(ind(ismem)))',[],@mode));

wholeCells = dendrites(wholeCellId);
dendrites = dendrites(setdiff(1:numel(dendrites),wholeCellId));
indBigDends = Agglo.isMaxBorderToBorderDistAbove(p, 5000, Superagglos.transformAggloNewOldRepr(dendrites));
[ myelinDend ] = connectEM.calculateSurfaceMyelinScore( dendrites, graph, borderMeta, heuristics ); % calculate myelin score for the dendrite class
    
end

