function splitChiasmaSub2(start_idx)
load('/tmpscratch/kboerg/chiasmarun/initial.mat');
load('/tmpscratch/kboerg/chiasmarun/initial_segids.mat');
thisEdgesNew = [];
thisEdgesCol = [];

for idx_top = start_idx : 500 :  length(idxtemp{which_col})
    load(['/tmpscratch/kboerg/chiasmarun/resultchiasma_' num2str(idx_top)]);
    idx_top
    [~,lookupnodes]=ismember(thisNodes,nodesScaled,'rows');
    assert(length(C) > 3 && sum(cellfun(@(idx)max(pdist2(thisNodes(idx, :), nodesScaled(i_here,:))) > 4000, C))>3);
    nrExits = length(C);
    if nrExits ~= 4
        continue
    end
    getIdsFromCC=@(x)skelsegids{1}(ismember(nodesScaled,thisNodes(x,:),'rows'));
    nonull=@(x)x(x~=0);
    whichones = find(cellfun(@(x)any(ismember(getIdsFromCC(x),nonull(ff.segIds{fflookup(idxtemp{which_col}(idx_top))}))),C));
    if length(whichones) ~= 2
        continue;
    end
    if isempty(thisEdgesCol)
        thisEdgesCol = lookupnodes(thisEdges);
    else
        thisEdgesCol=intersect(thisEdgesCol, lookupnodes(thisEdges), 'rows');
    end
    
    for idx = 1 : 2
        closest1 = pdist2(nodesScaled(lookupnodes(C{whichones(1)}),:), nodesScaled(temp{which_col}.output.ccCenterIdx(i_here),:));
        closest2 = pdist2(nodesScaled(lookupnodes(C{whichones(2)}),:), nodesScaled(temp{which_col}.output.ccCenterIdx(i_here),:));
        [~,X]=min(closest1);
        [~,Y]=min(closest2);
        
        thisEdgesNew = [thisEdgesNew; lookupnodes(C{whichones(1)}(X)),lookupnodes(C{whichones(2)}(Y))];
        temp{which_col}.output.nodes(lookupnodes(C{whichones(1)}(X)),:)
        temp{which_col}.output.nodes(lookupnodes(C{whichones(2)}(Y)),:)
        
        whichones = setdiff(1:4,whichones);
    end
end
save(['/tmpscratch/kboerg/chiasmarun2/resultchiasma2_' num2str(idx_top)],'thisEdgesNew','thisEdgesCol');