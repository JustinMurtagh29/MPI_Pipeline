function [newAgglos, summary] = ...
        splitChiasmataMulti(p, agglo, queries, varargin)
    % TODO(amotta):
    % • We still don't make use of the seed location / node of the
    %   tracings. In theory we only need to infer the end location.
    %
    % Written by
    %   Kevin Boergens <kevin.boergens@brain.mpg.de>
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    opts = struct;
    opts.outputFile = [];
    opts.dryRun = false;
    opts.exportNml = false;
    opts = Util.modifyStruct(opts, varargin{:});
    
    % configuration
    cutoffDistNm = 175;
    
    rng(0); % to make randperm further down reproducible
    
    % Group queries into chiasmata
   [~, chiasmata, queries.uniChiasmaId] = unique(queries.chiasmaId);
    chiasmata = queries(chiasmata, {'chiasmaId', 'centerNodeId'});
    chiasmaCount = size(chiasmata, 1);
    chiasmaTracings = accumarray( ...
        queries.uniChiasmaId, (1:size(queries, 1))', ...
        [], @(rows) {queries(rows, :)});
    
    % NOTE(amotta): Initialize key variables which are modified in the for
    % loop below. These variables collect all changes which need to be
    % applied to the original super-agglomerate.
    nodesToDelete = [];
    thisEdgesCol = agglo.edges;
    thisEdgesNew = zeros(0, 2);
    
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
        expectedNrExits = size(chiasmaTracings{chiIdx}, 1);
        
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
        
        % NOTE(amotta): Since we query all endings we expect the number of
        % tracings to match the number of found exits.
        nrExits = sum(isExit);
        assert(nrExits == expectedNrExits);
        
        summary.nrExits(chiIdx) = nrExits;
        summary.centerIdx(chiIdx) = centerIdx;
        summary.nrNonExits(chiIdx) = sum(~isExit);
        summary.nrExits(chiIdx) = nrExits;
        summary.tracings{chiIdx} = struct;
        summary.tracings{chiIdx}.nodes = ...
            chiasmaTracings{chiIdx}.flightNodes;
        summary.tracings{chiIdx}.overlaps = cell(nrExits, 1);
        
        % NOTE(amotta): Non-exit components are dropped (for now at least)
        nonExitNodeIds = thisNodeIds(cell2mat(C(~isExit)));
        C = C(isExit);
        
        % NOTE(amotta): Each entry of `groups` contains a list of exits
        % which were grouped together by virtue of chiasmata queries.
        groups = cell(1, 0);
        conns = zeros(0, 2);
        
        exitNodesScaled = chiasmaTracings{chiIdx}.seedPos;
        exitNodesScaled = bsxfun(@times, exitNodesScaled, p.voxelSize);
        
        for trIdx = 1:nrExits
            tr = chiasmaTracings{chiIdx}(trIdx, :);
            assert(tr.exitId == trIdx);
            
            % NOTE(amotta): To map a query to its chiasma exit we simply
            % walk from the tail of the flight path to its head. When we
            % first come into proximity of a (or multiple) exit nodes, we
            % stop. Proximity is defined by the `cutoffDistNm` distance
            % threshold. If there are multiple exit nodes within proximity,
            % then the closer one wins.
            
            % NOTE(amotta): A tracing may be empty (hence the reshape)
            tracingScaled = reshape(tr.flightNodes{1}, [], 3);
            tracingScaled = bsxfun(@times, tracingScaled, p.voxelSize);
            
            % NOTE(amotta): Make sure nodes are correctly sorted
            seedPosScaled = tr.seedPos .* p.voxelSize;
           [~, minIdx] = min(pdist2(seedPosScaled, tracingScaled));
            assert(minIdx == 1);
            
            % NOTE(amotta): Calculate distance from flight path to exit
            % nodes. But prevent overlap with seed by setting distance to
            % infinity.
            overlaps = pdist2(exitNodesScaled, tracingScaled);
            overlaps(tr.exitId, :) = inf;
            
            overlapNode = find(any(overlaps < cutoffDistNm, 1), 1, 'last');
           [~, overlaps] = min(overlaps(:, overlapNode));
           
            overlaps((end + 1):1) = 0;
            overlaps = cat(1, tr.exitId, overlaps);
            summary.tracings{chiIdx}.overlaps{trIdx} = overlaps;
        end
        
        if opts.dryRun; continue; end
        
        % NOTE(amotta): If we're unable to solve a chiasma, it is left
        % untouched and remains part of the super-agglomerate.
        if nrExits ~= sum(cellfun(@length, groups)) %unable to solve chiasma
            continue;
        end
        
        summary.solved(chiIdx) = true;
        
        % NOTE(amotta): Find index of node closest to chiasma center for
        % each component. These are the nodes from which the new edges will
        % be made.
        closest = nan(1, nrExits);
        for idx2 = 1 : nrExits
            [~,closest(idx2)] = min(pdist2( ...
                agglo.nodesScaled(thisNodeIds(C{idx2}), :), ...
                agglo.nodesScaled(centerIdx, :)));
            closest(idx2) = thisNodeIds(C{idx2}(closest(idx2)));
        end
        
        % NOTE(amotta): Translate connections to node IDs
        conns = sort(closest(conns), 2);
        
        % NOTE(amotta): Build skeleton where only core is cut out
        p.sphereRadiusOuter = Inf; % in nm
        [~, thisEdges, ~, thisNodeIds] = ...
             connectEM.detectChiasmataPruneToSphere( ...
             agglo.nodesScaled, agglo.edges, ...
             ones(size(agglo.edges, 1), 1), p, centerIdx);
        
        % NOTE(amotta): Non-exits are dropped from the agglomerate
        nodesToDelete = union(nodesToDelete, nonExitNodeIds);
        thisEdgesCol(any(ismember(thisEdgesCol, nonExitNodeIds), 2), :) = [];
        
        % NOTE(amotta): Build set of edges to keep. We start with the full
        % set and only keep the ones which never are part of a core.
        thisEdgesCol=intersect(thisEdgesCol, thisNodeIds(thisEdges), 'rows');
        nodesToDelete = union(nodesToDelete, setdiff( ...
            1:size(agglo.nodesScaled, 1), thisNodeIds));
        thisEdgesNew = [thisEdgesNew; conns];
    end
    
    % NOTE(amotta): See KMB's email from Wednesday, 27.09.2017.
    %
    % It's possible that the point-of-attachment for a solved chiasma query
    % is within another chiasma. As a result, the point-of-attachment is
    % pruned away and the new edge cannot be made. Is only one end of these
    % edges is pruned away, mark the remaining node as new open ending.
    prunedMask = any(ismember(thisEdgesNew, nodesToDelete), 2);
    newEndingNodes = setdiff(thisEdgesNew(prunedMask, :), nodesToDelete);
    thisEdgesNew(prunedMask, :) = [];
    
    % NOTE(amotta): Make sure that none of the nodes involved in edges is
    % being removed by one of the other chiasmata.
    assert(~any(ismember(thisEdgesCol(:), nodesToDelete)));
    assert(~any(ismember(thisEdgesNew(:), nodesToDelete)));
    
    %% build new super-agglomerates
    nodesToKeep = setdiff(1:size(agglo.nodesScaled,1), nodesToDelete);
    newMaxSegId = nodesToKeep(end);
    
    newEdges = cat(1, thisEdgesCol, thisEdgesNew);
    newNodeComps = Graph.findConnectedComponents(newEdges);
    newNodeComps = cat(1, newNodeComps, num2cell(reshape( ...
        setdiff(nodesToKeep, cell2mat(newNodeComps)), [], 1)));
    
    newNodeLUT = Agglo.buildLUT(newMaxSegId, newNodeComps);
    newEdgeComps = newNodeLUT(newEdges(:, 1));
    
    newAggloCount = numel(newNodeComps);
    newAgglos = struct;
    
    for curIdx = 1:newAggloCount
        curNodeIds = newNodeComps{curIdx};
        curNodes = agglo.nodes(curNodeIds, :);
        curEndings = find(ismember(curNodeIds, newEndingNodes));
        
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
