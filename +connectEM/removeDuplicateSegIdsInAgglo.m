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
        [~, idx] = unique(agglos(n).nodes(:,4),'first');  % this also works for NaNs (which are no duplicates
        idx = ~ismember(1:size(agglos(n).nodes,1),idx);  % get index to duplicates
        totalDupl = totalDupl + sum(idx);
        agglos(n).nodes(idx,:) = []; % remove overlapping nodes
        agglos(n).edges = unique(connectEM.changem(agglos(n).edges, 1:numel(idx)-cumsum(idx),1:numel(idx)),'rows');% correct edges accordingly
    end
end

