function splitChiasmataMulti(agglo, tasks, p, backup, aggloidx,outputFile)
    % Written by
    %   Kevin Boergens <kevin.boergens@brain.mpg.de>
    
    % TODO(amotta): Check where and how the `agglo` and `tasks` inputs are
    % generated. `agglo` must contain at least a `nodesScaled` field.
    % `tasks` is a structure array with `centeridx` containing the index of
    % the center node for chiasma `idx`.
    
    % TODO(amotta): Find out whether `nodesScaled` has a fourth column with
    % the segment ID.
    
    % TODO(amotta): Remove the `backup` input variable.
    % TODO(amotta): Find out whether this same code can now also be used
    % for the four-fold chiasmata.
    
    rng(aggloidx+1337); % to make randperm further down reproducible)
    
    % NOTE(amotta): Initialize key variables which are modified in the for
    % loop below. These variables collect all changes which need to be
    % applied to the original super-agglomerate.
    nodesToDelete = [];
    thisEdgesCol = agglo.edges;
    
    % TODO(amotta): Verify that sphere radius is the same one here as in
    % the processing of 4-fold chiasmata. That's what we want, I assume.
    p.voxelSize = [11.24 11.24 28];
    p.sphereRadiusInner = 1000; % in nm
    
    for idx = 1 : length(tasks)
        % NOTE(amotta): Restrict skeleton to components within shell
        p.sphereRadiusOuter = 10000; % in nm
        [thisNodes, thisEdges] = ... % cut down to sphere
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
        
        % NOTE(amotta): Make sure there are at least five connected
        % components. This must be true because this code was written to
        % handle >4 chiasmata.
        % TODO(amotta): Relax this to greater-equal if this same function
        % is also used to process the four-fold chiasmata.
        assert(length(C) > 4);
        nrExits = length(C);
        
        % TODO(amotta): Make sure each connected component has at least one
        % node that is more than 3 µm from the center node. Why?
        assert(all(cellfun( ...
            @(idx2) max(pdist2(thisNodes(idx2, :), ...
            agglo.nodesScaled(tasks(idx).centeridx,:))) > 3000, C)));
        
        % NOTE(amotta): Each entry of `groups` contains a list of exists
        % which were grouped together by virtue of chiasmata queries.
        groups = cell(1, 0);
        conns = zeros(0, 2);
        
        if nrExists == 4
            todoTracings = 1;
        else
            % TODO(amota): Find out why the upper limit is commented out
            todoTracings = 1; %: nrExits; % make list of tracings TEMPORARILY FOR 1st BATCH ONLY
        end
        
        % NOTE(amotta): Process flight paths until done. We're done when
        % one of the following conditions is true:
        %
        % * no more open exists
        % * no more flight paths to process
        % * exactly two open exits
        %
        % In the last case, the two endings are connected automatically.
        % This is the same behaviour as for the 4-fold chiasmata.
        while true
            nrOpenExits = nrExits - sum(cellfun(@length, groups));
            
            if nrOpenExits == 0
                break;
            elseif nrOpenExits == 2
                conns(end+1, :) = setdiff(1:nrExits, cell2mat(groups));
                groups(end+1) = {conns(end, :)};
                break;
            end
            
            % NOTE(amotta): `idx2pre` contains the list of open tracings
            randomOne = @(x)x(randperm(length(x),1));
            idx2pre = setdiff(todoTracings, cell2mat(groups));
            
            if ~isempty(idx2pre)
                idx2 = randomOne(idx2pre);
            elseif ~isempty(todoTracings)
                idx2 = randomOne(todoTracings);
            else
                break;
            end
            
            % NOTE(amotta): Mark tracing as processed
            todoTracings = setdiff(todoTracings, idx2);
            
            % NOTE(amotta): Find components touched by flight path
            getIdsFromCC=@(x)agglo.nodes(ismember(agglo.nodesScaled,thisNodes(x,:),'rows'),4);
            nonull = @(x)x(x~=0);
            whichones = find(cellfun(@(x)any(ismember(getIdsFromCC(x),nonull(tasks(idx).tracings(idx2).segids))),C))';
            
            if length(whichones) ~= 2
                % NOTE(amotta): Flight is invalid. Go to next...
                continue;
            end
            
            conns(end + 1, :) = whichones;
            groups = Graph.findConnectedComponents(conns);
        end
        
        % NOTE(amotta): If we're unable to solve a chiasma, it is left
        % untouched and remains part of the super-agglomerate.
        if nrExits ~= sum(cellfun(@length, groups)) %unable to solve chiasma
            continue;
        end
        
        % NOTE(amotta): Find index of node closest to chiasma center for
        % each component. These are the nodes from which the new edges will
        % be made.
        closest = nan(1, nrExists);
        for idx2 = 1 : nrExits
            [~,closest(idx2)] = min(pdist2( ...
                agglo.nodesScaled(lookupnodes(C{idx2}), :), ...
                agglo.nodesScaled(tasks(idx).centeridx, :)));
        end
        
        % NOTE(amotta): Build new edges
        thisEdgesNew = lookupnodes(closest(conns));
        thisEdgesNew = sort(thisEdgesNew, 2);
        
        % NOTE(amotta): Build skeleton where only core is cut out
        p.sphereRadiusOuter = Inf; % in nm
        [thisNodes, thisEdges] = ...
             connectEM.detectChiasmataPruneToSphere( ...
             agglo.nodesScaled, agglo.edges, ...
             ones(size(agglo.edges,1),1), p, tasks(idx).centeridx);
         
        C = Graph.findConnectedComponents(thisEdges);
        [~,lookupnodes]=ismember(thisNodes,agglo.nodesScaled,'rows');
        
        % NOTE(amotta): Build set of edges to keep. We start with the full
        % set and only keep the ones which never are part of a core.
        thisEdgesCol=intersect(thisEdgesCol, lookupnodes(thisEdges), 'rows');
        nodesToDelete = [nodesToDelete, setdiff(1:size(agglo.nodesScaled, 1), lookupnodes)];
    end
    
    % NOTE(amotta): Make sure that none of the nodes involved in edges is
    % being removed by one of the other chiasmata.
    assert(~any(ismember(thisEdgesCol(:), nodesToDelete)));
    assert(~any(ismember(thisEdgesNew(:), nodesToDelete)));
    
    agglo.nodes(nodesToDelete,:) = [];
    nodesToKeep = setdiff(1:size(agglo.nodesScaled,1), nodesToDelete);
    lookupNTK(nodesToKeep) = 1 : length(nodesToKeep);
    edgesToWrite = lookupNTK([thisEdgesCol;thisEdgesNew]);
    edgesToWrite = sort(edgesToWrite, 2);
    
    comments = cell(size(agglo.nodes,1), 1);
    %writeNml([outputFile(1:end-4) '.nml'], writeSkeletonFromNodesAndEdges({agglo.nodes}, {edgesToWrite}, {comments}, {'axon'}, {[0 0 1 1]}));
    
    newAgglo = struct;
    newAgglo.nodes = agglo.nodes;
    newAgglo.edges = edgesToWrite;
    Util.saveStruct(outputFile, newAgglo);
end
function dummy()
end