function [ newSuperAgglos ] = removeDuplicates( superAgglos )
% removes duplicate while making sure that edges are preserved

[~,keepIndices,uIndices] = arrayfun(@(x) unique(x.nodes(:,1:3),'rows'),superAgglos,'uni',0); %necessary to exclude 4th column because NaNs are always unique
[nodesToKeep, edgesToKeep] = arrayfun(@(x) deal(superAgglos(x).nodes(keepIndices{x},:),connectEM.changem(superAgglos(x).edges,uIndices{x},1:size(superAgglos(x).nodes,1))),1:numel(superAgglos),'uni',0);
newSuperAgglos = cell2struct([edgesToKeep';nodesToKeep'],{'edges','nodes'},1);

end

