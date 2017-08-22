function [ agglos, totalDupl] = removeDuplicateSegIdsInAgglo( agglos )
% this function removes duplicate segment ids in agglos (new or old
% representation)
%
% author: marcel.beining@brain.mpg.de

if ~isfield(agglos,'nodes') % old representation
    totalDupl = sum(cellfun(@numel,agglos));
    agglos = cellfun(@(x) unique(x),agglos,'uni',0);
    totalDupl = totalDupl - sum(cellfun(@numel,agglos));
else
    totalDupl = 0;
    for n = 1:numel(agglos)
        [~, idx,idx2] = unique(agglos(n).nodes(:,4),'first');  % this also works for NaNs (which are no duplicates
        agglos(n).nodes = agglos(n).nodes(idx,:); % keep only unique nodes
        agglos(n).edges = unique(connectEM.changem(agglos(n).edges, idx2,1:size(agglos(n).nodes,1)),'rows');% correct edges accordingly
        totalDupl = totalDupl + size(agglos(n).nodes,1) - numel(idx);
    end
end

