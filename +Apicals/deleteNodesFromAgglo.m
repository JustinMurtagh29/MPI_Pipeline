% delete nodes from agglo (new representation)
% in:   agglo (new rep) - agglomerate with edges and nodes
%       toDelete - logical vector, active for nodes to be deleted
% out:  newAgglo - updated edges and node list
% -------------------------------------------------------------------------
% Author: Matej Zecevic <matej.zecevic@brain.mpg.de>
function newAgglo = deleteNodesFromAgglo(agglo,toDelete)
tic;
% generate mask to take out all the nodes to be deleted
mask = ismember(agglo.edges, find(toDelete == 1));

% apply mask
newEdges = agglo.edges(all(mask,2),:);

% return new agglo
newAgglo.edges = newEdges;
newAgglo.nodes = agglo.nodes(~toDelete,:);

disp('finished execution of deleteNodesFromAgglo.');
toc; disp('---------------------------------------');

end

% agglosToWrite = [];
% for t=1:length(trees)
%    agglosToWrite = vertcat(agglosToWrite, {agglos{21}(vertcat(trees{t}))});
% end