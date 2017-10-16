function [newAgglos, summary] = ...
        splitChiasmataMulti(p, agglo, queries, varargin)
    % Written by
    %   Kevin Boergens <kevin.boergens@brain.mpg.de>
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    opts = struct;
    opts.outputFile = [];
    opts.dryRun = false;
    opts.exportNml = false;
    opts = Util.modifyStruct(opts, varargin{:});
    
    % configuration
    cutoffDistNm = 300; % only consider exits within X nm of flight path
    cutoffTracingNm = 1000; % ignore the first X nm of flight path
    
    % Group queries into chiasmata
   [~, chiasmata, queries.uniChiasmaId] = unique(queries.chiasmaId);
    chiasmata = queries(chiasmata, {'chiasmaId', 'centerNodeId'});
    chiasmaCount = size(chiasmata, 1);
    chiasmataTracings = accumarray( ...
        queries.uniChiasmaId, (1:size(queries, 1))', ...
        [], @(rows) {queries(rows, :)});
    
    % NOTE(amotta): Initialize key variables which are modified in the for
    % loop below. These variables collect all changes which need to be
    % applied to the original super-agglomerate.
    newEndings = [];
    nodesToDelete = [];
    edgesToKeep = agglo.edges;
    nodeCount = size(agglo.nodes, 1);
    
    nodesToAdd = zeros(0, 4);
    edgesToAdd = zeros(0, 2);
    
    summary = struct;
    summary.nrChiasmata = chiasmaCount;
    summary.centerIdx = nan(chiasmaCount, 1);
    summary.nrExits = nan(chiasmaCount, 1);
    summary.nrNonExits = nan(chiasmaCount, 1);
    summary.nrExits = nan(chiasmaCount, 1);
    summary.tracings = cell(chiasmaCount, 1);
    summary.solved = false(chiasmaCount, 1);
    
    p.voxelSize = p.raw.voxelSize;
    p.sphereRadiusInner = 1000; % in nm
    
    for chiIdx = 1:chiasmaCount
        centerIdx = chiasmata.centerNodeId(chiIdx);
        chiasmaTracings = chiasmataTracings{chiIdx};
        expectedNrExits = size(chiasmaTracings, 1);
        
        % NOTE(amotta): Restrict skeleton to components within shell
        p.sphereRadiusOuter = 10000; % in nm
       [thisNodes, thisEdges, ~, thisNodeIds] = ...
            connectEM.detectChiasmataPruneToSphere( ...
            agglo.nodesScaled, agglo.edges, ...
            ones(size(agglo.edges, 1), 1), p, centerIdx);
        
        C = Graph.findConnectedComponents(thisEdges);
        
        % NOTE(amotta): For a component to be considered an exit, it must
        % contain at least one node which is more than 3 µm from the
        % chiasma center.
        isExit = cellfun( ...
            @(idx2) max(pdist2(thisNodes(idx2, :), ...
            agglo.nodesScaled(centerIdx, :))) > 2000, C);
        
        % NOTE(amotta): Non-exit components are dropped (for now at least)
        nonExitNodeIds = thisNodeIds(cell2mat(C(~isExit)));
        nodesToDelete = union(nodesToDelete, nonExitNodeIds);
        edgesToKeep(any(ismember(edgesToKeep, nonExitNodeIds), 2), :) = [];
        
        % NOTE(amotta): Since we query all endings we expect the number of
        % tracings to match the number of found exits.
        C = C(isExit);
        nrExits = numel(C);
        assert(nrExits == expectedNrExits);
        
        summary.nrExits(chiIdx) = nrExits;
        summary.centerIdx(chiIdx) = centerIdx;
        summary.nrNonExits(chiIdx) = sum(~isExit);
        summary.nrExits(chiIdx) = nrExits;
        summary.tracings{chiIdx} = struct;
        summary.tracings{chiIdx}.nodes = ...
            chiasmaTracings.flightNodes;
        summary.tracings{chiIdx}.overlaps = cell(nrExits, 1);
        summary.tracings{chiIdx}.overlapNodes = nan(nrExits, 1);
        
        %%
        exitNodesScaled = chiasmaTracings.seedPos;
        exitNodesScaled = bsxfun(@times, exitNodesScaled, p.voxelSize);
        
        for trIdx = 1:nrExits
            tr = chiasmaTracings(trIdx, :);
            assert(tr.exitId == trIdx);
            
            % NOTE(amotta): To map a query to its chiasma exit we simply
            % search for the exit which comes closest to the flight path.
            % However, the first `cutoffTracingNm` nm of the flight path
            % are ignored. And so are exits which do not come closer than
            % `cutoffDistNm` to the flight path.
            
            % NOTE(amotta): A tracing may be empty (hence the reshape)
            trNodesScaled = reshape(tr.flightNodes{1}, [], 3);
            trNodesScaled = bsxfun(@times, trNodesScaled, p.voxelSize);
            
            % NOTE(amotta): Make sure nodes are correctly sorted
            trSeedScaled = tr.seedPos .* p.voxelSize;
           [~, minIdx] = min(pdist2(trSeedScaled, trNodesScaled));
            assert(isempty(trNodesScaled) || minIdx == 1);
            
            trDist = diff(trNodesScaled, 1, 1);
            trDist = sqrt(sum(trDist .* trDist, 2));
            trDist = cumsum(cat(1, 0, trDist)) < cutoffTracingNm;
            
            % NOTE(amotta): Calculate distance from flight path to exit
            % nodes. But prevent overlap with seed by setting distance to
            % infinity. Also ignore the first stretch of flight path.
            trProximity = pdist2(exitNodesScaled, trNodesScaled);
            trProximity(tr.exitId, :) = inf;
            trProximity(:, trDist) = inf;
            
           [trOverlapDist, overlaps] = min(trProximity(:));
           
            if trOverlapDist < cutoffDistNm
               [overlaps, overlapNode] =  ...
                   ind2sub(size(trProximity), overlaps);
            else
                overlapNode = 0;
                overlaps = 0;
            end
           
            overlaps = cat(1, tr.exitId, overlaps);
            summary.tracings{chiIdx}.overlaps{trIdx} = overlaps;
            summary.tracings{chiIdx}.overlapNodes(trIdx) = overlapNode;
        end
        
        if opts.dryRun; continue; end
        
        % NOTE(amotta): Build skeleton where only core is cut out
        p.sphereRadiusOuter = Inf; % in nm
       [~, thisEdges, ~, thisNodeIds] = ...
             connectEM.detectChiasmataPruneToSphere( ...
             agglo.nodesScaled, agglo.edges, ...
             ones(size(agglo.edges, 1), 1), p, centerIdx);
        
        % NOTE(amotta): Build set of edges to keep. We start with the full
        % set and only keep the ones which never are part of a core.
        edgesToKeep = intersect( ...
            edgesToKeep, thisNodeIds(thisEdges), 'rows');
        nodesToDelete = union( ...
            nodesToDelete, setdiff(1:nodeCount, thisNodeIds));
    end
    
    if opts.dryRun
        newAgglos = agglo;
        return;
    end
    
    %% decide how to solve chiasma
    % TODO(amotta): To some magic to determine which flight paths to
    % execute. This is determined in connectEM.evaluateChiasmataBatch
    % for now and then passed in here via tha `execute` variable of
    % `queries.`
        
    %% patch in flight paths
    for chiIdx = 1:chiasmaCount
        chiasmaTracings = chiasmataTracings{chiIdx};
        
        for trIdx = reshape(find(chiasmaTracings.execute), 1, [])
            tr = chiasmaTracings(trIdx, :);
            
            % segment ID `nan` on flight paths
            trNodes = tr.flightNodes{1};
            trNodes(:, 4) = nan;
            
            % cut off tail of flight path
            trLastNode = summary.tracings{chiIdx}.overlapNodes(trIdx);
            trNodes = trNodes(2:trLastNode, :);
            
            % collect new edges
            trNodeCount = size(trNodes, 1);
            trEdges = zeros((trNodeCount - 1) + 2, 2);
            trEdges((1 + 1):end, 1) = 1:trNodeCount;
            trEdges(1:(end - 1), 2) = 1:trNodeCount;
            
            trNodeOff = nodeCount + size(nodesToAdd, 1);
            trEdges = trEdges + trNodeOff;
            
            % find exit nodes to be connected
            trShortEdge = [ ...
                tr.exitNodeId, ...
                chiasmaTracings.exitNodeId(tr.overlaps{1}(end))];
            trEdges([1, end]) = trShortEdge;
            
            % It's possible that one of the two exit nodes got already
            % pruned away by cutting out one of the earlier chiasmata. For
            % more info, see KMB's email from Wednesday, 27.09.2017
            trShortEdge = setdiff(trShortEdge, nodesToDelete);
            
            if numel(trShortEdge) == 2
                % path in flight path
                nodesToAdd = cat(1, nodesToAdd, trNodes);
                edgesToAdd = cat(1, edgesToAdd, trEdges);
            else
                % mark new endings, but discard flight path
                newEndings = union(newEndings, trShortEdge);
            end
        end
    end
    
    % NOTE(amotta): Make sure that none of the nodes involved in edges is
    % being removed by one of the other chiasmata.
    assert(~any(ismember(edgesToKeep(:), nodesToDelete)));
    assert(~any(ismember(edgesToAdd(:), nodesToDelete)));
    
    %% build new super-agglomerates
    % build new nodes
    newNodes = cat(1, agglo.nodes, nodesToAdd);
    nodesToKeep = setdiff(1:size(newNodes, 1), nodesToDelete);
    newNodes = newNodes(nodesToKeep, :);
    
    % renumber node IDs in `newEndings`
   [~, newEndings] = ismember(newEndings, nodesToKeep);
    newEndings(~newEndings) = [];
    
    % renumber node IDs in edges
    newEdges = cat(1, edgesToKeep, edgesToAdd);
   [~, newEdges] = ismember(newEdges, nodesToKeep);
    assert(all(newEdges(:)));
   
    newNodeComps = Graph.findConnectedComponents(newEdges);
    newNodeComps = cat(1, newNodeComps, num2cell(reshape( ...
        setdiff(size(newNodes, 1), cell2mat(newNodeComps)), [], 1)));
    
    newNodeLUT = Agglo.buildLUT(size(newNodes, 1), newNodeComps);
    newEdgeComps = newNodeLUT(newEdges(:, 1));
    
    newAggloCount = numel(newNodeComps);
    newAgglos = struct;
    
    for curIdx = 1:newAggloCount
        curNodeIds = newNodeComps{curIdx};
        curNodes = newNodes(curNodeIds, :);
        curEndings = find(ismember(curNodeIds, newEndings));
        
        curEdges = newEdges(newEdgeComps == curIdx, :);
        [~, curEdges] = ismember(curEdges, curNodeIds);
        assert(all(curEdges(:)));
        
        newAgglos(curIdx).nodes = curNodes;
        newAgglos(curIdx).edges = curEdges;
        newAgglos(curIdx).endings = curEndings;
    end
    
    newAgglos = reshape(newAgglos, [], 1);
    
    % Save results
    if ~isempty(opts.outputFile)
        Util.save( ...
            strcat(opts.outputFile, '.mat'), ...
            agglo, queries, newAgglos, summary);
    end
    
    % NML for debuggin
    if opts.exportNml
        skel = skeleton();
        
        comments = repmat({''}, size(agglo.nodes, 1), 1);
        comments(cat(1, tasks.centeridx)) = {'Chiasma'};
        
        skel = skel.addTree( ...
            'Original', agglo.nodes, agglo.edges, [], [], comments);
        skel = skel.addBranchpoint(cat(1, tasks.centeridx));
        
        for curIdx = 1:newAggloCount
            skel = skel.addTree( ...
                sprintf('Component %d', curIdx), ...
                newAgglos(curIdx).nodes, newAgglos(curIdx).edges);
        end
        
        skel = Skeleton.setParams4Pipeline(skel, p);
        skel.write(strcat(opts.outputFile, '.nml'));
        clear comments skel;
    end
end
