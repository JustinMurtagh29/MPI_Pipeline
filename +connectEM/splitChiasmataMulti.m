function [newAgglos, summary] = ...
        splitChiasmataMulti(p, agglo, tasks, outputFile)
    % Written by
    %   Kevin Boergens <kevin.boergens@brain.mpg.de>
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % configuration
    doExportNml = false;
    
    rng(0); % to make randperm further down reproducible
    
    % NOTE(amotta): Initialize key variables which are modified in the for
    % loop below. These variables collect all changes which need to be
    % applied to the original super-agglomerate.
    nodesToDelete = [];
    thisEdgesCol = agglo.edges;
    thisEdgesNew = zeros(0, 2);
    
    summary = struct;
    summary.nrChiasmata = numel(tasks);
    summary.nrExits = nan(numel(tasks), 1);
    summary.nrNonExits = nan(numel(tasks), 1);
    summary.nrTracings = nan(numel(tasks), 1);
    summary.tracings = cell(numel(tasks), 1);
    summary.solved = false(numel(tasks), 1);
    
    p.voxelSize = p.raw.voxelSize;
    p.sphereRadiusInner = 1000; % in nm
    
    for idx = 1 : length(tasks)
        % NOTE(amotta): Restrict skeleton to components within shell
        p.sphereRadiusOuter = 10000; % in nm
        [thisNodes, thisEdges, ~, thisNodeIds] = ...
            connectEM.detectChiasmataPruneToSphere( ...
            agglo.nodesScaled, agglo.edges, ...
            ones(size(agglo.edges,1),1), p, tasks(idx).centeridx);
        
        % make sure we have the same conn components as for the detection 
        % so that our indexing of tasks(idx).tracings(idx2) works
        % i.e. so that tracings(idx2) starts at C{idx2}
        %assert(isequal(thisNodes, backup.thisNodes)); NO CHECKING RIGHT NOW
        %assert(isequal(thisEdges, backup.thisEdges)); NO CHECKING RIGHT NOW
        
        C = Graph.findConnectedComponents(thisEdges);
        
        % NOTE(amotta): For a component to be considered an exit, it must
        % contain at least one node which is more than 3 µm from the
        % chiasma center.
        isExit = cellfun( ...
            @(idx2) max(pdist2(thisNodes(idx2, :), ...
            agglo.nodesScaled(tasks(idx).centeridx,:))) > 3000, C);
        
        % NOTE(amotta): Make sure there are at least four connected
        % components. This must be true because this code was written to
        % handle >=4 chiasmata.
        assert(sum(isExit) >= 4);
        nrExits = sum(isExit);
        nrTracings = numel(tasks(idx).tracings);
        
        summary.nrExits(idx) = nrExits;
        summary.nrNonExits(idx) = sum(~isExit);
        summary.nrTracings(idx) = nrTracings;
        summary.tracings{idx} = struct;
        summary.tracings{idx}.processed = zeros(nrTracings, 1);
        summary.tracings{idx}.overlaps = cell(nrTracings, 1);
        
        if nrExits == 4
            todoTracings = 1;
        else
            todoTracings = 1:nrExits;
        end
        
        % NOTE(amotta): Non-exit components are dropped (for now at least)
        nonExitNodeIds = thisNodeIds(cell2mat(C(~isExit)));
        C = C(isExit);
        
        % NOTE(amotta): Each entry of `groups` contains a list of exits
        % which were grouped together by virtue of chiasmata queries.
        groups = cell(1, 0);
        conns = zeros(0, 2);
        tracingIdx = 1;
        
        % NOTE(amotta): Process flight paths until done. We're done when
        % one of the following conditions is true:
        %
        % * no more open exits
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
            
            summary.tracings{idx}.processed(idx2) = tracingIdx;
            summary.tracings{idx}.overlaps{idx2} = whichones;
            
            if length(whichones) ~= 2
                % NOTE(amotta): Flight is invalid. Go to next...
                continue;
            end
            
            conns(end + 1, :) = whichones;
            groups = Graph.findConnectedComponents(conns);
            tracingIdx = tracingIdx + 1;
        end
        
        % NOTE(amotta): If we're unable to solve a chiasma, it is left
        % untouched and remains part of the super-agglomerate.
        if nrExits ~= sum(cellfun(@length, groups)) %unable to solve chiasma
            continue;
        end
        
        summary.solved(idx) = true;
        
        % NOTE(amotta): Find index of node closest to chiasma center for
        % each component. These are the nodes from which the new edges will
        % be made.
        closest = nan(1, nrExits);
        for idx2 = 1 : nrExits
            [~,closest(idx2)] = min(pdist2( ...
                agglo.nodesScaled(thisNodeIds(C{idx2}), :), ...
                agglo.nodesScaled(tasks(idx).centeridx, :)));
            closest(idx2) = thisNodeIds(C{idx2}(closest(idx2)));
        end
        
        % NOTE(amotta): Translate connections to node IDs
        conns = sort(closest(conns), 2);
        
        % NOTE(amotta): Build skeleton where only core is cut out
        p.sphereRadiusOuter = Inf; % in nm
        [~, thisEdges, ~, thisNodeIds] = ...
             connectEM.detectChiasmataPruneToSphere( ...
             agglo.nodesScaled, agglo.edges, ...
             ones(size(agglo.edges,1),1), p, tasks(idx).centeridx);
        
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
    
    % NOTE(amotta): The rest is only executed if `outputFile` is set
    if ~exist('outputFile', 'var') || isempty(outputFile); return; end
    Util.save(strcat(outputFile, '.mat'), agglo, tasks, newAgglos, summary);
    
    %% for debugging
    if doExportNml
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
        skel.write(strcat(outputFile, '.nml'));
        clear comments skel;
    end
end
