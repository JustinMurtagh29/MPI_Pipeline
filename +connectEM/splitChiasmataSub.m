function splitChiasmataSub(idx_start)
load('/tmpscratch/kboerg/chiasmarun/initial.mat');
for idx_top = idx_start : 2000: length(idxtemp{which_col})
    i_here = find(cellfun(@(x)any(ismember(x,ff.startNode{fflookup(idxtemp{1}(idx_top))},'rows')),temp{1}.output.position));
    lost = length(i_here) -1;
    i_here=i_here(1);
    %currently we can't deal if two tracings start at the same start node. We should keep track of the task ids for that
    [thisNodes, thisEdges, thisProb] = ...
        connectEM.detectChiasmataPruneToSphere( ...
        nodesScaled, temp{which_col}.output.edges, ...
        ones(size(temp{which_col}.output.edges,1),1), p, temp{which_col}.output.ccCenterIdx(i_here));
    C = Graph.findConnectedComponents(thisEdges);
    save(['/tmpscratch/kboerg/chiasmarun/resultchiasma_' num2str(idx_top)], 'C', 'thisNodes','thisEdges','i_here','lost');
end