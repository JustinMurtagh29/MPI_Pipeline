function splitChiasmataSub(idx_start)
    load('/tmpscratch/kboerg/chiasmarunAugust/preSplitRun');
    for idx = idx_start:50:length(idx_pre)
        superagglo_definition = load(['/tmpscratch/kboerg/visX23/visX19_' num2str(floor(idx_pre(idx)/100)) '/visX19_' num2str(idx_pre(idx)) '/result.mat']);
        copyfile(['/tmpscratch/kboerg/visX23/visX19_' num2str(floor(idx_pre(idx)/100)) '/visX19_' num2str(idx_pre(idx)) '/skelDetected.nml'],[outputFolder 'skelDetectedOld' num2str(idx) '.nml']);
        
        

        % get skels' ids
        skelsegids = connectEM.lookupSkelGT(p, struct('nodes',{{superagglo_definition.output.nodes}}));
        nodesScaled = bsxfun(@times,superagglo_definition.output.nodes,p.voxelSize);
        selected = lookup_tasks(lookup_tasks(:,1)==idx_pre(idx)&idx_pre_pre',2);
        thisEdgesNew = [];
        nodesToDelete = [];
        thisEdgesCol = superagglo_definition.output.edges;
        for idx2 = 1 : length(selected);
            p.sphereRadiusOuter = 10000; % in nm
            p.sphereRadiusInner = 1000; % in nm
            p.voxelSize = [11.24 11.24 28];
            [thisNodes, thisEdges, thisProb] = ...
                connectEM.detectChiasmataPruneToSphere( ...
                nodesScaled, superagglo_definition.output.edges, ...
                ones(size(superagglo_definition.output.edges,1),1), p, superagglo_definition.output.ccCenterIdx(selected(idx2)));
            C = Graph.findConnectedComponents(thisEdges);
            [~,lookupnodes]=ismember(thisNodes,nodesScaled,'rows');
            assert(length(C) > 3 && sum(cellfun(@(idx)max(pdist2(thisNodes(idx, :), nodesScaled(superagglo_definition.output.ccCenterIdx(selected(idx2)),:))) > 3000, C))>3);
            nrExits = length(C);
            ff_idx = find(cellfun(@length,strfind(ff.filenames,M(ismember(lookup_tasks,[idx_pre(idx),selected(idx2),1],'rows')))));
            copyfile(ff.filenames{ff_idx},[outputFolder 'query' num2str(idx) '_' num2str(idx2) '.nml']);
            if nrExits ~= 4
                continue
            end
            getIdsFromCC=@(x)skelsegids{1}(ismember(nodesScaled,thisNodes(x,:),'rows'));
            nonull=@(x)x(x~=0);
            whichones = find(cellfun(@(x)any(ismember(getIdsFromCC(x),nonull(ff.segIds{ff_idx}))),C));
            if length(whichones) ~=2
                continue
            end
            
            
            for idx3 = 1 : 2
                closest1 = pdist2(nodesScaled(lookupnodes(C{whichones(1)}),:), nodesScaled(superagglo_definition.output.ccCenterIdx(selected(idx2)),:));
                closest2 = pdist2(nodesScaled(lookupnodes(C{whichones(2)}),:), nodesScaled(superagglo_definition.output.ccCenterIdx(selected(idx2)),:));
                [~,X]=min(closest1);
                [~,Y]=min(closest2);
                
                thisEdgesNew = [thisEdgesNew; lookupnodes(C{whichones(1)}(X)),lookupnodes(C{whichones(2)}(Y))];
                superagglo_definition.output.nodes(lookupnodes(C{whichones(1)}(X)),:)
                superagglo_definition.output.nodes(lookupnodes(C{whichones(2)}(Y)),:)
                
                whichones = setdiff(1:4,whichones);
            end
            p.sphereRadiusOuter = Inf; % in nm
            p.sphereRadiusInner = 1000; % in nm
            p.voxelSize = [11.24 11.24 28];
            [thisNodes, thisEdges, thisProb] = ...
                connectEM.detectChiasmataPruneToSphere( ...
                nodesScaled, superagglo_definition.output.edges, ...
                ones(size(superagglo_definition.output.edges,1),1), p, superagglo_definition.output.ccCenterIdx(selected(idx2)));
            C = Graph.findConnectedComponents(thisEdges);
            [~,lookupnodes]=ismember(thisNodes,nodesScaled,'rows');
            thisEdgesCol=intersect(thisEdgesCol, lookupnodes(thisEdges), 'rows');
            nodesToDelete = [nodesToDelete, setdiff(1:size(nodesScaled, 1), lookupnodes)];
            
        end
        superagglo_definition.output.nodes(nodesToDelete,:) = [];
        nodesToKeep = setdiff(1:size(nodesScaled,1), nodesToDelete);
        lookupNTK(nodesToKeep) = 1 : length(nodesToKeep);
        edgesToWrite = lookupNTK([thisEdgesCol;thisEdgesNew]);
        comments = cell(size(superagglo_definition.output.nodes,1), 1);
        writeNml([outputFolder 'skelSplit' num2str(idx) '.nml'], writeSkeletonFromNodesAndEdges({superagglo_definition.output.nodes}, {edgesToWrite}, {comments}, {'axon'}, {[0 0 1 1]}));
        
    end    
end
function dummy()
end