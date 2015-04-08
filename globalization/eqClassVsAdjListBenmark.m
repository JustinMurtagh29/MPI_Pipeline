function eqClassVsAdjListBenmark(filenameInput)

    load(filenameInput);
   
    display('Evaluating CC calculation methods: (debug & optimize findConnectedComponents): ');
    tic;
    [a,b] = findConnectedComponents(graph1.edges);
    [c,d] = findConnectedComponents(graph2.edges);
    toc;
    
%     display('First approach (): ');
%     % Join those two mappings into one 
%     bothMappings = [graph1.eqClass; graph2.eqClass];
%     clear components allComponents;
%     % Initalize empty adjaceny list of equivalence classes
%     adjList = [];
%     lengthBothMappings = length(bothMappings);
%     for i=1:lengthBothMappings
%         temp = find(cellfun(@(x,y)any(ismember(x,y)), repmat(bothMappings(i),lengthBothMappings-i,1), bothMappings(i+1:end)))+i;
%         adjList = [adjList; i.*ones(length(temp),1) temp];
%     end
%     bothMappingsIdx = findConnectedComponents(adjList, unique(adjList));
%     oneMapping = cell(length(bothMappingsIdx),1);
%     for i=1:length(bothMappingsIdx)
%         oneMapping{i} = bothMappings{bothMappingsIdx{i}};
%     end



end
