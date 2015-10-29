function oneMappingForProblemGeneration(p, upperT)
    % Mainly for testing GP based joining of supervoxel graph, maybe should become class later
    % Or at least put more of this stuff into more modular functions
    % WORK on this function seized after noticing that two components are two orders of magnitude LARGER than the others
    % TO DO: Find out why!

    if ~exist([p.saveFolder 'graph.mat'], 'file')
        tic;
        % Collect all edges of local supervoxel graph
        graph.edges = NaN(1e8,2); % 1e8 heuristic for 07x2, turns out actually ~7.5e7
        graph.prob = NaN(1e8,1);
        graph.cubeLI = NaN(1e8,1); % store linear indices as well (second column for correspondences)
        idx = 1;
        % IMPORTANT TO DO: Remove index shifts here in next run, just quick fix for now
        for i=2:size(p.local,1)-1
            for j=1:size(p.local,2)
                for k=1:size(p.local,3)
                    load(p.local(i,j,k).edgeFile);
                    load(p.local(i,j,k).probFile); 
                    % Put everything in one structure
                    nrElements = length(prob);
                    graph.edges(idx:idx+nrElements-1,:) = edges;
                    graph.prob(idx:idx+nrElements-1) = prob;
                    graph.cubeLI(idx:idx+nrElements-1) = repmat(sub2ind(size(p.local),i,j,k), [nrElements 1]);
                    idx = idx+nrElements;
                end
            end
        end

        % Collect all correspondences between local supervoxel graph
        files = dir([p.correspondence.saveFolder, '*global.mat']);
        for i=1:length(files)
            % IMPORTANT TO DO: same here, change to using back all indices!
            if ~(strcmp(files(i).name(1:2),'01') || strcmp(files(i).name(1:2), '15') || ...
                    strcmp(files(i).name(7:8),'01') || strcmp(files(i).name(7:8), '15'))
                load([p.correspondence.saveFolder files(i).name]);
                nrElements = length(corrGlobal1);
                graph.edges(idx:idx+nrElements-1,:) = [corrGlobal1 corrGlobal2];
                graph.prob(idx:idx+nrElements-1) = ones(length(corrGlobal1),1);
                idx = idx+nrElements;
            end
        end
        % Drop part of arrays preallocated but not assigned
        graph.edges(idx:end,:) = [];
        graph.prob(idx:end) = [];
        graph.cubeLI(idx:end) = [];
        save([p.saveFolder 'graph.mat'], 'graph');
        toc;
    else
        display('NOTE: Graph loaded from HD. Rename for recalculation');
        load([p.saveFolder 'graph.mat']);
    end

    % Dont do this right now, big components, perculation etc. (no idea) :)
    if ~exist([p.saveFolder 'graphNew.mat'], 'file');
        tic;
        % Find maximum ID before joining to determine new (joined) supervoxel IDs
        load([p.local(end,end,end).saveFolder 'localToGlobalSegId.mat']);
        maxGlobalId = max(globalIds);
        % Upper threshold
        idxJoined = graph.prob > upperT;
        % Some bookkepping
        edgesJoined = graph.edges(idxJoined,:);
        %graph.probJoined = graph.prob(idxJoined);
        %graph.edgesRemaining = graph.edges(~idxJoined,:);
        %graph.probRemaining = graph.prob(~idxJoined);
        [graph.ccEdgesJoined.equivalenceClasses, graph.ccEdgesJoined.objectClassLabels] = findConnectedComponents(edgesJoined);
        % Loop for relabeling edges and calculating probabilities accordingly
        if 0 % Do not globalize graph for now (components membership is checked in gameProblemTesting)
            for i=1:length(graph.ccEdgesJoined.equivalenceClasses)
                display(num2str(i));
                tic;
                % Find all remaining edges pointing into current (i-th) connected component (not within, see all in 2nd & 3rd line)
                logicalFindEdges = ismember(graph.edgesRemaining, graph.ccEdgesJoined.equivalenceClasses{i});
                idxIntoComponent = any(logicalFindEdges,2) & ~all(logicalFindEdges,2);
                idxWithinComponent = all(logicalFindEdges,2);
                % Store edges and probabilities into component before removing to calculate what to (re)insert
                graph.ccEdgesJoined.edgesInto{i} = graph.edgesRemaining(idxIntoComponent,:);
                graph.ccEdgesJoined.probInto{i} = graph.probRemaining(idxIntoComponent);
                graph.ccEdgesJoined.edgesWithin{i} = graph.edgesRemaining(idxWithinComponent,:);
                graph.ccEdgesJoined.probWithin{i} = graph.probRemaining(idxWithinComponent);           
                % Remove all (into and within) of those from remaining edges and probabilities
                graph.edgesRemaining(idxIntoComponent|idxWithinComponent,:) = [];
                graph.probRemaining(idxIntoComponent|idxWithinComponent) = [];
                % Renumber all edges into component to new supervoxel ID 
                edgesNew = graph.ccEdgesJoined.edgesInto{i};
                edgesNew(logicalFindEdges(idxIntoComponent,:)) = maxGlobalId + i;
                % Make sure first edge column is always smaller
                edgesNew = sort(edgesNew,2);
                % Remove double entries / redundancies (thereby generating new edges into component after mapping)   
                [graph.ccEdgesJoined.edgesNew{i}, idxNew, idxOld] = unique(edgesNew, 'rows');
                % Calculate new probabilities
                idxMatrix = sparse(1:length(idxOld), idxOld, 1);
                probMatrix = bsxfun(@times, idxMatrix, graph.ccEdgesJoined.probInto{i});
                probNew = zeros([size(probMatrix,2) 1]);
                for j=1:length(probNew)
                    probNew(j) = 1 - prod(1 - full(nonzeros(probMatrix(:,j))));
                end
                graph.ccEdgesJoined.probNew{i} = probNew;
                % Reinsert into graph
                graph.edgesRemaining = [graph.edgesRemaining; graph.ccEdgesJoined.edgesNew{i}];
                graph.probRemaining = [graph.probRemaining; graph.ccEdgesJoined.probNew{i}];
                toc;
            end
        end
        % Save for future use 
        save([p.saveFolder 'graphNew.mat'], 'graph', '-v7.3');
        toc;
    else
        display('Graph after GP and correspodence based joining of supervoxel loaded from disk, NOT recalculated')
        load([p.saveFolder 'graphNew.mat']);
    end

end

