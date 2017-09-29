function splitChiasmataSub(idx_start)
    load('/tmpscratch/kboerg/chiasmarunAugust/preSplitRun');
    
    % NOTE(amotta): Iterate over agglomerates
    for idx = idx_start:50:length(idx_pre)
        % NOTE(amotta): `idx_pre(idx)` is the index of the current agglo
        superagglo_definition = load(['/tmpscratch/kboerg/visX23/visX19_' num2str(floor(idx_pre(idx)/100)) '/visX19_' num2str(idx_pre(idx)) '/result.mat']);
        copyfile(['/tmpscratch/kboerg/visX23/visX19_' num2str(floor(idx_pre(idx)/100)) '/visX19_' num2str(idx_pre(idx)) '/skelDetected.nml'],[outputFolder 'skelDetectedOld' num2str(idx) '.nml']);
        
        % NOTE(amotta): Look up segment IDs at nodes of agglomerate.
        %
        % There significantly faster ways to do this. Also, this is
        % reconstructing discarded information. But code should be right.
        skelsegids = connectEM.lookupSkelGT(p, struct('nodes',{{superagglo_definition.output.nodes}}));
        nodesScaled = bsxfun(@times,superagglo_definition.output.nodes,p.voxelSize);
        
        % NOTE(amotta): Find the indices of the answered chiasmata!
        selected = lookup_tasks(lookup_tasks(:,1)==idx_pre(idx)&idx_pre_pre',2); %#ok
        
        % NOTE(amotta): Initalize the three variables which will be used to
        % construct the resulting agglomerates.
        thisEdgesNew = [];
        nodesToDelete = [];
        thisEdgesCol = superagglo_definition.output.edges;
        assert(issorted(thisEdgesCol, 2));
        
        % NOTE(amotta): Iterate over chiasmata
        for idx2 = 1 : length(selected)
            p.sphereRadiusOuter = 10000; % in nm
            p.sphereRadiusInner = 1000; % in nm
            p.voxelSize = [11.24 11.24 28];
            
            % NOTE(amotta): Removes the connected component within the
            % sphere centered on the current node.
            [thisNodes, thisEdges] = ...
                connectEM.detectChiasmataPruneToSphere( ...
                nodesScaled, superagglo_definition.output.edges, ...
                ones(size(superagglo_definition.output.edges,1),1), p, ...
                superagglo_definition.output.ccCenterIdx(selected(idx2)));
            
            % NOTE(amotta): Map node subset back to whole agglomerate
            [~,lookupnodes]=ismember(thisNodes,nodesScaled,'rows');
            
            C = Graph.findConnectedComponents(thisEdges);
            
            % NOTE(amotta): By definition a chiasma has at least four
            % exits, such that the sphere pruning yields at least four
            % different connected components.
            nrExits = length(C);
            assert(length(C) > 3);
            
            % NOTE(amotta): Make sure that in each connected component
            % there is at least one node with more than 3 µm distance from
            % the current center node.
            %
            % Not yet sure why this is required...
            
            % TODO(amotta): Shouldn't this be `== nrExists` instead of
            % `> 3`? Anyway, it's correct for four-fold chiasmata.
            assert(sum(cellfun(@(idx) max(pdist2(thisNodes(idx, :), nodesScaled(superagglo_definition.output.ccCenterIdx(selected(idx2)),:))) > 3000, C))>3);
            
            % NOTE(amotta): Find indices of focused flight tasks for
            % current chiasma.
            ff_idx = ismember(lookup_tasks,[idx_pre(idx),selected(idx2),1],'rows');
            ff_idx = find(cellfun(@length,strfind(ff.filenames, M(ff_idx))));
            
            copyfile(ff.filenames{ff_idx},[outputFolder 'query' num2str(idx) '_' num2str(idx2) '.nml']);
            
            % NOTE(amotta): Only process four-fold chiasmata for now
            if nrExits ~= 4
                continue
            end
            
            % NOTE(amotta): Find overlap between the 4+ connected
            % components and the flight path. `whichones` contains the
            % indices of the two components to be reconnected.
            getIdsFromCC=@(x)skelsegids{1}(ismember(nodesScaled,thisNodes(x,:),'rows'));
            nonull=@(x)x(x~=0);
            whichones = find(cellfun(@(x)any(ismember(getIdsFromCC(x),nonull(ff.segIds{ff_idx}))),C));
            
            % NOTE(amotta): Only apply chiasma splitting if exactly two
            % components were touched by the flight path. This is correct
            % for the case of a four-fold chiasma which consists of two
            % merged neurites.
            if length(whichones) ~=2
                continue
            end
            
            % For a pair of agglomerates to-connect
            % * find the nodes closest to the chiasma center
            % * add a new edge between them
            for idx3 = 1 : 2
                closest1 = pdist2(nodesScaled(lookupnodes(C{whichones(1)}),:), nodesScaled(superagglo_definition.output.ccCenterIdx(selected(idx2)),:));
                closest2 = pdist2(nodesScaled(lookupnodes(C{whichones(2)}),:), nodesScaled(superagglo_definition.output.ccCenterIdx(selected(idx2)),:));
                [~,X]=min(closest1);
                [~,Y]=min(closest2);
                
                thisEdgesNew = [thisEdgesNew; lookupnodes(C{whichones(1)}(X)),lookupnodes(C{whichones(2)}(Y))]; %#ok
                
                % NOTE(amotta): Show nodes to be connected
                superagglo_definition.output.nodes(lookupnodes(C{whichones(1)}(X)),:)
                superagglo_definition.output.nodes(lookupnodes(C{whichones(2)}(Y)),:)
                
                % NOTE(amotta): Switch to other pair
                % TODO(amotta): This only works for four-fold chiasmata
                whichones = setdiff(1:4,whichones);
            end
            
            p.sphereRadiusOuter = Inf; % in nm
            p.sphereRadiusInner = 1000; % in nm
            p.voxelSize = [11.24 11.24 28];
            
            [thisNodes, thisEdges] = ...
                connectEM.detectChiasmataPruneToSphere( ...
                nodesScaled, superagglo_definition.output.edges, ...
                ones(size(superagglo_definition.output.edges,1),1), p, ...
                superagglo_definition.output.ccCenterIdx(selected(idx2)));
            assert(issorted(thisEdges, 2));
            
            [~,lookupnodes]=ismember(thisNodes,nodesScaled,'rows');
            assert(issorted(lookupnodes(thisEdges), 2));
            
            thisEdgesCol=intersect(thisEdgesCol, lookupnodes(thisEdges), 'rows');
            nodesToDelete = [nodesToDelete, setdiff(1:size(nodesScaled, 1), lookupnodes)];
            
        end
        
        % NOTE(amotta): Results from above are:
        % * A list of edges to keep
        % * A list of edges to add
        % * A list of segments to remove
        %
        % Both of these use indices to the **original** agglomerate. After
        % removing the nodes, the edge-to-add must thus be renumbered.
        superagglo_definition.output.nodes(nodesToDelete,:) = [];
        
        nodesToKeep = setdiff(1:size(nodesScaled,1), nodesToDelete);
        lookupNTK = zeros(1, nodesToKeep(end));
        lookupNTK(nodesToKeep) = 1 : length(nodesToKeep);
        
        edgesToWrite = lookupNTK([thisEdgesCol;thisEdgesNew]);
        edgesToWrite = sort(edgesToWrite, 2);
        
        % Save NML file
        comments = cell(size(superagglo_definition.output.nodes,1), 1);
        writeNml([outputFolder 'skelSplit' num2str(idx) '.nml'], writeSkeletonFromNodesAndEdges({superagglo_definition.output.nodes}, {edgesToWrite}, {comments}, {'axon'}, {[0 0 1 1]}));
        
    end    
end
function dummy()
end