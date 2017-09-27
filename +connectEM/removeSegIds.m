function [ agglo ] = removeSegIds(agglo,segIds)
%removes segment ids from an agglo (new representation) without checking
%for split agglos

for n = 1:numel(agglo)
    idxToDelete = find(ismember(agglo(n).nodes(:,4),segIds))';
    omap = (1:size(agglo(n).nodes,1))';
    newmap = omap-sum(bsxfun(@ge,omap,idxToDelete),2);
    newmap(idxToDelete) = 0;
    agglo(n).edges = connectEM.changem(agglo(n).edges,newmap,omap);
    if all(agglo(n).edges(:)==0) % to avoid that goddamn bug in matlab which makes 0x1 instead of 0x2 edge vector
        agglo(n).edges = zeros(0,2);
    else
        agglo(n).edges(any(agglo(n).edges==0,2),:) = [];
    end
    agglo(n).nodes(idxToDelete,:) = [];
end

end

