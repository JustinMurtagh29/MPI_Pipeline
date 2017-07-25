function [nodes2, edges2] = makeSkelForChiasmataDetectionSub(currentEC, axons, edgesGTall, resultCol)
last = @(x)x{end};
fifth = @(x)x{5}
getTask = @(x)fifth(strsplit(last(strsplit(x,'/')),'_'));

nodes = [];
edges = [];
lookup1 = [];
lookup2 = [];
% create nodes and edges within agglos
for idx = 1 : length(currentEC)
    idx
    lookup=sparse(ones(1,length(axons{currentEC(idx)})), axons{currentEC(idx)}, 1:length(axons{currentEC(idx)}));
    
    edges=[edges;full(lookup(edgesGTall(all(ismember(edgesGTall,axons{currentEC(idx)}),2),:)))+size(nodes,1)];
    nodes=[nodes;segmentMeta.point(:,axons{currentEC(idx)})'];
    lookup1 = [lookup1; repmat(currentEC(idx), length(axons{currentEC(idx)}), 1)];
    lookup2 = [lookup2; axons{currentEC(idx)}];
end
edges2=edges;
nodes2=nodes;
usedTasks = {};
% create nodes and edges for queries
for runidx = 1 : 2
    for idx = 1 : length(resultCol{runidx}.startAgglo)
        if mod(idx,100)==0
            disp(idx);
        end
        if ~isempty(resultCol{runidx}.startAgglo{idx}) && resultCol{runidx}.startAgglo{idx}<=length(axons) && ismember(resultCol{runidx}.startAgglo{idx},currentEC)
            if size(resultCol{runidx}.ff.nodes{idx},1)<2 ||any(ismember(getTask(resultCol{runidx}.ff.filenames{idx}),usedTasks))
                continue;
            end
            % make sure each task is only used once (discarding quality control redundancy
            usedTasks{end+1}=getTask(resultCol{runidx}.ff.filenames{idx});
            % find nodes that connect to agglos
            tempids =[resultCol{runidx}.ff.segIds{idx},resultCol{runidx}.ff.neighbours{idx}];
            %somehow we lost the node order for the query, here reconstructed with MSP
            Tree = graphminspantree(sparse(squareform(pdist(resultCol{runidx}.ff.nodes{idx}))));
            [X,Y]=find(Tree);
            [~, Locb] = ismember(tempids(:),lookup2);
            [I, ~] = ind2sub(size(tempids),find(Locb));
            edges2=[edges2;[X,Y]+size(nodes2,1); I+size(nodes2,1), Locb(Locb>0)];
            nodes2=[nodes2;resultCol{runidx}.ff.nodes{idx}];
            assert(size(nodes2,1)>=max(edges2(:)));
        end
    end
end
