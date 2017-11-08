function superagglos = get50RandDend(superagglos)
% returns the 50 pseudorandom agglos for any agglo state

load(fullfile('/gaba/u/mberning/results/pipeline/20170217_ROI/aggloState','randDendIds50.mat'),'segIdsForRandDend')

[aggloLUT,aggloSegIds] = Superagglos.buildLUT(agglos);

[~,ind] = ismember(segIdsForRandDend,aggloSegIds);
randDends = NaN(numel(segIdsForRandDend),1);
randDends(ind~=0) = aggloLUT(aggloSegIds(ind(ind~=0)));

superagglos = superagglos(randDends);