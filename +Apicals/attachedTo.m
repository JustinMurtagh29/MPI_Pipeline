% attachedTo 
% for a set of segments (agglo) get segments to which spine heads
% attached to (designed for L4)
% -------------------------------------------------------------------------
% Inputs:
% agglo - dendrite agglo (designed for components, subsets)
% segmentsReached - list of spine head -> segment
% shCounts, k (optional) - spine head counts & index of agglo
% Outputs:
% shAttached - for a specific agglo the segments where sh attached to
% -------------------------------------------------------------------------
% Author: Matej Zecevic <matej.zecevic@brain.mpg.de>
function [ shAttached ] = attachedTo( agglo, segmentsReached, k, shCounts )

% get segments spine heads attached to
shAttached = [];
count = 0;
for s=1:length(agglo) 
    ind = find(agglo(s) == segmentsReached);
    if length(ind) >= 1
        for g=1:length(ind)
            shAttached = vertcat(shAttached, segmentsReached(ind(g)));
        end
    elseif ~isempty(ind)
        shAttached = vertcat(shAttached, segmentsReached(ind));
    end
end

% check if result is consistent
% number of spine heads attached == total times segments reached
if exist('k','var') && exist('shCounts','var')
    assert(sum(ismember(shAttached, agglo)) == shCounts(k));
    disp('completed assertion successfully.');
end


end

