function [ newclass1, equivalencesClass1 ] = executeEdges(class1,edges,segmentMeta)
% executes the edges between agglos of class1, and also adds
% single segments of edges

allEdges = cell2mat(arrayfun(@(x) reshape(class1(x).nodes(class1(x).edges,4),[],2),1:numel(class1),'uni',0)');
if all(all(ismember(edges,allEdges),2))
    newclass1 = class1;
    return
else
    allEdges = cat(1,allEdges, edges);
    
    [equivalencesClass1, labelObj] = Graph.findConnectedComponents(allEdges,0,1);
    
    newnodes = cellfun(@(x) [segmentMeta.point(x,:)',x],equivalencesClass1);
    [~, newedges] = cellfun(@(x) ismember(allEdges(labelObj(allEdges(:,1))==x,:),newnodes{x}(:,4)),1:numel(equivalencesClass1),'uni',0);
    newclass1 = cell2struct([newedges';newnodes'],{'edges','nodes'},1);
    
end

