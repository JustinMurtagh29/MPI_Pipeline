function oneMappingForProblemGeneration(p, upperT)
% Mainly for testing GP based joining of supervoxel graph, maybe should become class later

if ~exist([p.saveFolder 'graph.mat'], 'file')
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
else
    display('NOTE: Graph loaded from HD. Rename for recalculation');
    load([p.saveFolder 'graph.mat']);
end

if ~exist([p.saveFolder 'graphNew.mat'], 'file');
    % Upper threshold
    idxJoined = graph.prob > upperT;
    graph.edgesJoined = graph.edges(idxJoined,:);
    graph.probJoined = graph.prob(idxJoined);
    % Find connected components in joined edges 
    [graph.ccEdgesJoined.equivalenceClasses, graph.ccEdgesJoined.objectClassLabels] = findConnectedComponents(graph.edgesJoined);
    % Recalculate edge and probability list for joined edges
    % Initalize with edges not joined
    graph.edgesRemaining = graph.edges(~idxJoined,:);
    graph.probRemaining = graph.prob(~idxJoined);
    tempEdges = permute(graph.edgesRemaining, [3 1 2 ]);
    % Loop for relabeling edges and calculating probabilities accordinglyy
    for i=1:length(graph.ccEdgesJoined.equivalenceClasses)
        % Find all edges pointing into current (i-th) connected component
        tempFindEdges = bsxfun(@eq, tempEdges, graph.ccEdgesJoined.equivalenceClasses{i});
        edges = graph.edgesRemaining(squeeze(sum(sum(tempFindEdges,1),3)),:); 
        graph.edgesRemaining = a;
        graph.probRemaining = b;
    end
    % Save for future use 
    save([p.saveFolder 'graphNew.mat'], 'graph');
else
    display('Graph after GP and correspodence based joining of supervoxel loaded from disk, NOT recalculated')
    load([p.saveFolder 'graphNew.mat']);
end

end

