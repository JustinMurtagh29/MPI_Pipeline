function oneMappingForProblemGeneration(p, upperT)
% Mainly for testing GP based joining of supervoxel graph, maybe should become class later

if ~exist([p.saveFolder 'graph.mat'], 'file')
    tic;
    % Collect all edges of local supervoxel graph
    graph.edges = NaN(1e8,2); % 1e8 heuristic for 07x2, turns out actually ~7.5e7
    graph.prob = NaN(1e8,1);
    idx = 1;
    for i=1:size(p.local,1)
        for j=1:size(p.local,2)
            for k=1:size(p.local,3)
                load(p.local(i,j,k).edgeFile);
                load(p.local(i,j,k).probFile); 
                % Put everything in one structure
                nrElements = length(prob);
                graph.edges(idx:idx+nrElements-1,:) = edges;
                graph.prob(idx:idx+nrElements-1) = prob;
                idx = idx+nrElements;
            end
        end
    end

    % Collect all correspondences between local supervoxel graph
    files = dir([p.correspondence.saveFolder, '*global.mat']);
    for i=1:length(files)
        load([p.correspondence.saveFolder files(i).name]);
        nrElements = length(corrGlobal1);
        graph.edges(idx:idx+nrElements-1,:) = [corrGlobal1 corrGlobal2];
        graph.prob(idx:idx+nrElements-1) = ones(length(corrGlobal1),1);
        idx = idx+nrElements;
    end
    % Drop part of arrays preallocated but not assigned
    graph.edges(idx:end,:) = [];
    graph.prob(idx:end) = [];
    save([p.saveFolder 'graph.mat'], 'graph');
    toc;
else
    display('NOTE: Graph loaded from HD. Rename for recalculation');
    load([p.saveFolder 'graph.mat']);
end

if ~exist([p.saveFolder 'graphNew.mat'], 'file');
    tic;
    % Find maximum ID before joining to determine new (joined) supervoxel IDs
    load([p.local(end,end,end).saveFolder 'localToGlobalSegId.mat']);
    maxGlobalId = max(globalIds);
    % Upper threshold
    idxJoined = graph.prob > upperT;
    graph.edgesJoined = graph.edges(idxJoined,:);
    graph.probJoined = graph.prob(idxJoined);
    graph.edgesRemaining = graph.edges(~idxJoined,:);
    graph.probRemaining = graph.prob(~idxJoined);
    % Find connected components in joined edges 
    [graph.ccEdgesJoined.equivalenceClasses, graph.ccEdgesJoined.objectClassLabels] = findConnectedComponents(graph.edgesJoined);
    % Loop for relabeling edges and calculating probabilities accordingly
    for i=1:length(graph.ccEdgesJoined.equivalenceClasses)
        % Find all edges pointing into current (i-th) connected component (not within, see all in last line)
        logicalFindEdges = ismember(graph.edgesRemaining, graph.ccEdgesJoined.equivalenceClasses{i});
        idxToRemove = any(logicalFindEdges,2) & ~all(logicalFindEdges,2);
        % Store edges and probabilities before removing to calculate what to (re)insert
        edgesToRemove = graph.edgesRemaining(idxToRemove,:);
        probToRemove = graph.probRemaining(idxToRemove);
        % Remove all of those from remaining edges and probabilities
        graph.edgesRemaining(idxToRemove,:) = [];
        graph.probRemaining(idxToRemove) = [];
        % Renumber all edges to new supervoxel ID 
        edgesToRemove(logicalFindEdges(idxToRemove,:)) = maxGlobalId + i;
        % Make sure first edge column is always smaller
        edgesToRemove = sort(edgesToRemove,2);
        % Remove double entries / redundancies (thereby generating new edges after thresholod)   
        [edgesNew, ~, idxOld] = unique(edgesToRemove, 'rows');
        % Calculate new probabilities
        probNew = zeros(size(edgesToRemove,1),1);
        for j=1:size(edgesToRemove,1)
            probNew(j) = 1 - prod(1 - probToRemove(idxOld == j));
        end 
        % Reinsert into graph
        graph.edgesRemaining = [graph.edgesRemaining; edgesNew];
        graph.probRemaining = [graph.probRemaining; probNew];
    end
    % Save for future use 
    save([p.saveFolder 'graphNew.mat'], 'graph');
    toc;
else
    display('Graph after GP and correspodence based joining of supervoxel loaded from disk, NOT recalculated')
    load([p.saveFolder 'graphNew.mat']);
end

end

