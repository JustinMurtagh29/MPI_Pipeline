function splitChiasmataMulti(agglo, tasks, p, backup, aggloidx,outputFile)
    rng(aggloidx+1337); % to make randperm further down reproducible)
    thisEdgesNew = [];
    nodesToDelete = [];
    p.voxelSize = [11.24 11.24 28];
    p.sphereRadiusInner = 1000; % in nm
    
    thisEdgesCol = agglo.edges;
    for idx = 1 : length(tasks)
        p.sphereRadiusOuter = 10000; % in nm
        [thisNodes, thisEdges, thisProb] = ... % cut down to sphere
            connectEM.detectChiasmataPruneToSphere( ...
            agglo.nodesScaled, agglo.edges, ...
            ones(size(agglo.edges,1),1), p, tasks(idx).centeridx);
        % make sure we have the same conn components as for the detection 
        % so that our indexing of tasks(idx).tracings(idx2) works
        % i.e. so that tracings(idx2) starts at C{idx2}
        %assert(isequal(thisNodes, backup.thisNodes)); NO CHECKING RIGHT NOW
        %assert(isequal(thisEdges, backup.thisEdges)); NO CHECKING RIGHT NOW
        
        C = Graph.findConnectedComponents(thisEdges);
        [~,lookupnodes]=ismember(thisNodes,agglo.nodesScaled,'rows');
        assert(length(C) > 3 && sum(cellfun(@(idx2)max(pdist2(thisNodes(idx2, :), agglo.nodesScaled(tasks(idx).centeridx,:))) > 3000, C))>3);
        nrExits = length(C);
        groups = {};
        conns = [];
        if length(C) == 4
            todoTracings = 1;
        else
            todoTracings = 1 %: nrExits; % make list of tracings TEMPORARILY FOR 1st BATCH ONLY
        end
        while true
            if nrExits == sum(cellfun(@length, groups)) % none open
                break
            end
            if nrExits == sum(cellfun(@length, groups)) + 2 % only two open
                groups(end+1) = {setdiff(1:nrExits, cell2mat(groups))};
                break
            end
             
            randomOne = @(x)x(randperm(length(x),1));
            idx2pre = setdiff(todoTracings, cell2mat(groups));
            if isempty(idx2pre) % all obvious used
                if isempty(todoTracings) % all used
                    break
                else
                    idx2 = randomOne(todoTracings)
                end
            else
                idx2 = randomOne(idx2pre);
            end
            todoTracings = setdiff(todoTracings, idx2);
            
            getIdsFromCC=@(x)agglo.nodes((ismember(agglo.nodesScaled,thisNodes(x,:),'rows')),4);
            nonull = @(x)x(x~=0);
            whichones = find(cellfun(@(x)any(ismember(getIdsFromCC(x),nonull(tasks(idx).tracings(idx2).segids))),C))';
            if length(whichones)==2 % connection found
                conns = [conns; whichones];
            end
            if ~isempty(conns)
                 groups = Graph.findConnectedComponents(conns); %make new groups
            end     
        end
        if nrExits ~= sum(cellfun(@length, groups)) %unable to solve chiasma
            continue;
        end
        for idx2 = 1 : nrExits
            [~,closest(idx2)]=min(pdist2(agglo.nodesScaled(lookupnodes(C{idx2}),:), agglo.nodesScaled(tasks(idx).centeridx,:)));
        end
        for idx2 = 1 : length(groups) % add all the edges
            for idx3 = 1 : length(groups{idx2})
                for idx4 = (idx3 + 1) : length(groups{idx2});
                    thisEdgesNew = [thisEdgesNew; closest(groups{idx2}(idx3)),closest(groups{idx2}(idx4))];
                end
            end
        end
        p.sphereRadiusOuter = Inf; % in nm
        [thisNodes, thisEdges, thisProb] = ... % find nodes outside of inner sphere
             connectEM.detectChiasmataPruneToSphere( ...
             agglo.nodesScaled, agglo.edges, ...
             ones(size(agglo.edges,1),1), p, tasks(idx).centeridx);
        C = Graph.findConnectedComponents(thisEdges);
        [~,lookupnodes]=ismember(thisNodes,agglo.nodesScaled,'rows');
        thisEdgesCol=intersect(thisEdgesCol, lookupnodes(thisEdges), 'rows');
        nodesToDelete = [nodesToDelete, setdiff(1:size(agglo.nodesScaled, 1), lookupnodes)];
    end
    assert(~any(any(ismember([thisEdgesCol;thisEdgesNew],nodesToDelete))));
    agglo.nodes(nodesToDelete,:) = [];
    nodesToKeep = setdiff(1:size(agglo.nodesScaled,1), nodesToDelete);
    lookupNTK(nodesToKeep) = 1 : length(nodesToKeep);
    edgesToWrite = lookupNTK([thisEdgesCol;thisEdgesNew]);
    comments = cell(size(agglo.nodes,1), 1);
    %writeNml([outputFile(1:end-4) '.nml'], writeSkeletonFromNodesAndEdges({agglo.nodes}, {edgesToWrite}, {comments}, {'axon'}, {[0 0 1 1]}));
    nodes2 = agglo.nodes;
    edges2 = edgesToWrite;

    save(outputFile,'nodes2','edges2');
end
function dummy()
end