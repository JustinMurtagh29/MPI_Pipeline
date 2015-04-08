function oneMappingForProblemGeneration(p)

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
    % Load cube border correspondences and all local cube components based on GP
    load([p.saveFolder 'globalCorrespondences.mat']);
    load([p.saveFolder 'joinedComponents.mat']);
    % Join those two mappings into one 
    bothMappings = [components; allComponents];
    clear components allComponents;
    % Initalize empty adjaceny list of equivalence classes
    adjList = [];
    lengthBothMappings = length(bothMappings);
    for i=1:lengthBothMappings
        temp = find(cellfun(@(x,y)any(ismember(x,y)), repmat(bothMappings(i),lengthBothMappings-i,1), bothMappings(i+1:end)))+i;
        adjList = [adjList; i.*ones(length(temp),1) temp];
    end
    bothMappingsIdx = findConnectedComponents(adjList, unique(adjList));
    oneMapping = cell(length(bothMappingsIdx),1);
    for i=1:length(bothMappingsIdx)
        oneMapping{i} = bothMappings{bothMappingsIdx{i}};
    end
    save([p.saveFolder 'oneMapping.mat']);
else
    warning('Joining quotient sets of both mappings loaded from disk, NOT recalculated')
    load([p.saveFolder 'bothMappings.mat']);
end
% Iterate over cubes and apply thresholds to GP prediction
for i=1:size(p.local,1)
    for j=1:size(p.local,2)
        for k=1:size(p.local,3)
            load(p.local(i,j,k).edgeFile);
            load(p.local(i,j,k).probFile);
            % Drop edges over upper threshold (already CC's detected in joinedComponents.mat, variable componentsNew)
            edges = edges(prob < upperT);
            prob = prob(prob < upperT);
            % Drop edges below lower threshold
            edges = edges(prob > lowerT);
            prob = prob(prob > lowerT);j
            % Renumber edges according joinedComponents.mat, 
            for i=1:length(componentsNew)
                edges == componentsNew{i}(j);
            end
            % Put everything in one structure
            graph(i,j,k) = [edges prob];
        end
    end
end


% Save global correspodences


end
