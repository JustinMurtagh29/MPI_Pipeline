function [newAgglos, summary] = ...
        splitChiasmataMulti(p, chiParam, agglo, queries, varargin)
    % Written by
    %   Kevin Boergens <kevin.boergens@brain.mpg.de>
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    opts = struct;
    opts.outputFile = [];
    opts.dryRun = false;
    opts.exportNml = false;
    opts.partialAnswers = false;
    
    opts = Util.modifyStruct(opts, varargin{:});
    
    % configuration
    cutoffDistNm = 300; % only consider exits within X nm of flight path
    cutoffTracingNm = 1000; % ignore the first X nm of flight path
    cutoffExitNm = chiParam.minNodeDist; % discard exits shorter than X nm
    
    % override defaults
    if isfield(chiParam, 'split')
        if isfield(chiParam.split, 'cutoffDistNm') ...
                && ~isempty(chiParam.split.cutoffDistNm)
            cutoffDistNm = chiParam.split.cutoffDistNm;
        end
        
        if isfield(chiParam.split, 'cutoffTracingNm') ...
                && ~isempty(chiParam.split.cutoffTracingNm)
            cutoffTracingNm = chiParam.split.cutoffTracingNm;
        end
    end
    
    axonId = unique(queries.axonId);
    assert(isscalar(axonId));
    
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
    nodesToDelete = cell(chiasmaCount, 1);
    edgesToDelete = cell(chiasmaCount, 1);
    nodeCount = size(agglo.nodes, 1);
    
    nodesToAdd = zeros(0, 4);
    edgesToAdd = zeros(0, 2);
    
    summary = struct;
    summary.axonId = axonId;
    summary.nrChiasmata = chiasmaCount;
    summary.chiasmaId = nan(chiasmaCount, 1);
    summary.centerIdx = nan(chiasmaCount, 1);
    summary.nrExits = nan(chiasmaCount, 1);
    summary.nrNonExits = nan(chiasmaCount, 1);
    summary.solved = false(chiasmaCount, 1);
    summary.tracings = cell(chiasmaCount, 1);
    
    p.sphereRadiusInner = chiParam.sphereRadiusInner; % in nm
    
    for chiIdx = 1:chiasmaCount
        chiasmaId = chiasmata.chiasmaId(chiIdx);
        centerIdx = chiasmata.centerNodeId(chiIdx);
        chiTracings = chiasmataTracings{chiIdx};
        
        chiTracingIds = find(chiTracings.flightId);
        chiNrTracings = numel(chiTracingIds);
        
        % NOTE(amotta): Restrict skeleton to components within shell
        p.sphereRadiusOuter = chiParam.sphereRadiusOuter; % in nm
       [thisNodes, thisEdges, thisNodeIds] = ...
            connectEM.detectChiasmataPruneToSphere( ...
                p, agglo.nodesScaled, agglo.edges, centerIdx);
        
        C = Graph.findConnectedComponents( ...
            thisEdges, false, false, numel(thisNodeIds));
        
        % NOTE(amotta): For a component to be considered an exit, it must
        % contain at least one node which is more than 3 µm from the
        % chiasma center.
        isExit = cellfun( ...
            @(idx2) max(pdist2(thisNodes(idx2, :), ...
            agglo.nodesScaled(centerIdx, :))) > cutoffExitNm, C);
        
        % NOTE(amotta): Non-exit components are dropped (for now at least)
        nonExitNodeIds = thisNodeIds(cell2mat(C(~isExit)));
        
        C = C(isExit);
        nrExits = numel(C);
        
        expectedNrExits = size(chiTracings, 1);
        assert(nrExits == expectedNrExits);
        
        summary.chiasmaId(chiIdx) = chiasmaId;
        summary.centerIdx(chiIdx) = centerIdx;
        summary.nrExits(chiIdx) = nrExits;
        summary.nrNonExits(chiIdx) = sum(~isExit);
        summary.tracings{chiIdx} = struct;
        summary.tracings{chiIdx}.taskIds = ...
            chiTracings.taskId(chiTracingIds);
        summary.tracings{chiIdx}.nodes = ...
            chiTracings.flightNodes(chiTracingIds);
        summary.tracings{chiIdx}.overlaps = cell(chiNrTracings, 1);
        summary.tracings{chiIdx}.overlapNodes = nan(chiNrTracings, 1);
        summary.tracings{chiIdx}.execute = false(chiNrTracings, 1);
        
        %%
        exitNodesScaled = chiTracings.seedPos;
        exitNodesScaled = bsxfun(@times, exitNodesScaled, p.raw.voxelSize);
        
        for trIdx = 1:chiNrTracings
            trId = chiTracingIds(trIdx);
            tr = chiTracings(trId, :);
            
            % NOTE(amotta): To map a query to its chiasma exit we simply
            % search for the exit which comes closest to the flight path.
            % However, the first `cutoffTracingNm` nm of the flight path
            % are ignored. And so are exits which do not come closer than
            % `cutoffDistNm` to the flight path.
            
            % NOTE(amotta): A tracing may be empty (hence the reshape)
            trNodesScaled = reshape(tr.flightNodes{1}, [], 3);
            trNodesScaled = bsxfun(@times, trNodesScaled, p.raw.voxelSize);
            
            % NOTE(amotta): Make sure nodes are correctly sorted
            trSeedScaled = tr.seedPos .* p.raw.voxelSize;
           [~, minIdx] = min(pdist2(trSeedScaled, trNodesScaled));
           
            if ~(isempty(trNodesScaled) || minIdx == 1)
                warning('First node does not coincide with seed');
            end
            
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
       [~, thisEdges, thisNodeIds] = ...
             connectEM.detectChiasmataPruneToSphere( ...
                p, agglo.nodesScaled, agglo.edges, centerIdx);
        
        % NOTE(amotta): Find node and edges to delete. There are two
        % contributors to this list:
        % * cutting out of the inner sphere
        % * short branchlets which leaves the pruned sphere
        nodesToDelete{chiIdx} = cat( ...
            1, nonExitNodeIds(:), ...
            setdiff((1:nodeCount)', thisNodeIds));
        edgesToDelete{chiIdx} = find( ...
            any(ismember(agglo.edges, nonExitNodeIds), 2) ...
         | ~ismember(agglo.edges, thisNodeIds(thisEdges), 'rows'));
    end
    
    %% decide which chiasmata to split
    summary = connectEM.splitChiasmataMultiLogic(summary);
    
    if opts.dryRun
        newAgglos = agglo;
        return;
    end
    
    %% assemble final list of nodes and edges to remove
    nodesToDelete = cell2mat(nodesToDelete(summary.split));
    edgesToDelete = cell2mat(edgesToDelete(summary.split));
    
    %% patch in flight paths
    for chiIdx = 1:chiasmaCount
        chiTracings = chiasmataTracings{chiIdx};
        
        % check if chiasma is scheduled for splitting
        if ~summary.split(chiIdx); continue; end
        
        chiSummary = summary.tracings{chiIdx};
        chiExecuteIds = find(chiSummary.execute);
        
        % mark 1-components as new open endings
        chiOpenExitNodes = find(accumarray(chiSummary.lut, 1) == 1);
        chiOpenExitNodes = ismember(chiSummary.lut, chiOpenExitNodes);
        chiOpenExitNodes = chiTracings.exitNodeId(chiOpenExitNodes);
        newEndings = union(newEndings, chiOpenExitNodes);
        
        % determine short edges
        chiShortEdges = chiSummary.overlaps(chiExecuteIds);
        chiShortEdges = reshape(chiShortEdges, 1, []);
        chiShortEdges = transpose(cell2mat(chiShortEdges));
        chiShortEdges = chiTracings.exitNodeId(chiShortEdges);

        if any(ismember(chiShortEdges(:), nodesToDelete))
            % We cannot attach the flight paths and thus not split the
            % chiasma if any of the needed exit nodes got pruned away.
            % Instead, the remaining nodes are marked as new open endings.
            newEndings = union( ...
                newEndings, setdiff(chiShortEdges, nodesToDelete));
            continue;
        end
        
        for trIdx = 1:numel(chiExecuteIds)
            trId = chiExecuteIds(trIdx);
            tr = chiTracings(trId, :);
            
            % segment ID `nan` on flight paths
            trNodes = tr.flightNodes{1};
            trNodes(:, 4) = nan;
            
            % cut off tail of flight path
            trLastNode = chiSummary.overlapNodes(trId);
            trNodes = trNodes(2:trLastNode, :);
            
            % collect new edges
            trNodeCount = size(trNodes, 1);
            trEdges = zeros((trNodeCount - 1) + 2, 2);
            trEdges((1 + 1):end, 1) = 1:trNodeCount;
            trEdges(1:(end - 1), 2) = 1:trNodeCount;
            
            trNodeOff = nodeCount + size(nodesToAdd, 1);
            trEdges = trEdges + trNodeOff;
            
            % find exit nodes to be connected
            trEdges([1, end]) = chiShortEdges(trIdx, :);
            
            nodesToAdd = cat(1, nodesToAdd, trNodes);
            edgesToAdd = cat(1, edgesToAdd, trEdges);
            
            % mark flight path as executed
            summary.tracings{chiIdx}.execute(trId) = true;
        end
        
        % mark chiasma as solved
        summary.solved(chiIdx) = true;
    end
    
    % mark center nodes of solved chiasmata as solved
    agglo.solvedChiasma(summary.centerIdx(summary.solved)) = true;
    
    % NOTE(amotta): Make sure that none of the nodes involved in edges is
    % being removed by one of the other chiasmata.
    assert(~any(ismember(edgesToAdd(:), nodesToDelete)));
    
    %% build new super-agglomerates
    % build new nodes
    newNodes = cat(1, agglo.nodes, nodesToAdd);
    nodesToKeep = setdiff(1:size(newNodes, 1), nodesToDelete);
    newNodes = newNodes(nodesToKeep, :);
    
    % build new solved chiasma
    newSolvedChiasma = cat(1, ...
        agglo.solvedChiasma, false(size(nodesToAdd, 1), 1));
    newSolvedChiasma = newSolvedChiasma(nodesToKeep);
    
    % renumber node IDs in `newEndings`
   [~, newEndings] = ismember(newEndings, nodesToKeep);
    newEndings(~newEndings) = [];
    
    % renumber node IDs in edges
    newEdges = cat(1, agglo.edges, edgesToAdd);
    newEdges(edgesToDelete, :) = [];
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
        curSolvedChiasma = newSolvedChiasma(curNodeIds);
        curEndings = find(ismember(curNodeIds, newEndings));
        
        curEdges = newEdges(newEdgeComps == curIdx, :);
        [~, curEdges] = ismember(curEdges, curNodeIds);
        assert(all(curEdges(:)));
        
        
        newAgglos(curIdx).nodes = curNodes;
        newAgglos(curIdx).edges = curEdges;
        newAgglos(curIdx).endings = curEndings;
        newAgglos(curIdx).solvedChiasma = curSolvedChiasma;
    end
    
    newAgglos = reshape(newAgglos, [], 1);
    
    % Save results
    if ~isempty(opts.outputFile)
        Util.save( ...
            strcat(opts.outputFile, '_axonId_', num2str(axonId), '.mat'), ...
            agglo, queries, newAgglos, summary);
    end
    
    % NML for debugging
    if opts.exportNml
        skel = skeleton();
        
        comments = repmat({''}, size(agglo.nodes, 1), 1);
        comments(cat(1, summary.centerIdx)) = {'Chiasma'};
        
        skel = skel.addTree( ...
            'Original', agglo.nodes, agglo.edges, [], [], comments);
        skel = skel.addBranchpoint(cat(1, summary.centerIdx));
        
        for curIdx = 1:newAggloCount
            skel = skel.addTree( ...
                sprintf('Component %d', curIdx), ...
                newAgglos(curIdx).nodes, newAgglos(curIdx).edges);
        end
        
        skel = Skeleton.setParams4Pipeline(skel, p);
        skel.write(strcat(opts.outputFile, '_axonId_', num2str(axonId), '.nml'));
        clear comments skel;
    end
end
