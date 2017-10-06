function skels = newAxonGTOverlapsToSkel( gtSkels, sagglos, ov, minOv )
%NEWAXONGTOVERLAPSTOSKEL Write the gt axons and the overlapping
% agglomerates to skeletons.
% INPUT gtSkels: [10x1] cell
%           The new_axon_gt skeletons.
%           (see also connectEM.eval.getNewAxonGT)
%       sagglos: struct
%           The axon agglos in the superagglo format.
%       ov: [10x1] cell of [Nx2] int
%           The linear indices of superagglos that overlap with the
%           corresponding gtSkels in the first column and the number of
%           overlapping segments in the second column.
%           (see also connectEM.eval.getNewAxonGTAggloOverlap)
%       minOV: (Optional) double
%           Minimal number of overlapping segments to consider the
%           skeleton/agglo pair as overlapping.
%           (Default: 1)
% OUTPUT skels: [10x1] cell
%           Cell array of skeletons containing the ground truth skeletons
%           and the overlapping agglos.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('minOv', 'var') || isempty(minOv)
    minOv = 1;
end

skels = gtSkels;
for i = 1:length(skels)
    Util.log('Processing axon %d.', i);
    toAddAgglos = ov{i}(ov{i}(:,2) >= minOv, 1);
    skels{i} = L4.Agglo.superAgglo2Skel(sagglos(toAddAgglos), skels{i});
    skels{i}.names(2:end) = arrayfun(@(x)sprintf('Agglo_%d', x), ...
        toAddAgglos, 'uni', 0);
end

end

