function eqClassVsAdjListBenmark(filenameInput)

    load(filenameInput);

    % First check connected component implementation vs. matlab-version vs. bioinformatics version
    display('Evaluate CC calculation methods: ');
    tic;
    graph1.adjMatrix = sparse( graph1.edges(:,1), graph1.edges(:,2), 1, max(graph1.edges(:)), max(graph1.edges(:)) );
    graph1.adjMatrix = graph1.adjMatrix + graph1.adjMatrix';
    [S C] = conncomp(graph1.adjMatrix);
    graph1.segmentsToEqClass = C;
    graph1.eqClassToSegments = cell(S,1);
    bsxfun(unique(C));
    toc;

    display('First approach (): ');
    % Join those two mappings into one 
    bothMappings = [graph1.eqClass; graph2.eqClass];
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



end

function [S,C] = conncomp(G)
    [p,q,r] = dmperm(G'+speye(size(G)));
    S = numel(r)-1;
    C = cumsum(full(sparse(1,r(1:end-1),1,1,size(G,1))));
    C(p) = C;
end

