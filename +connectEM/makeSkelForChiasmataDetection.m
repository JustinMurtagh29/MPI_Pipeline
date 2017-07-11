borderGT = []
for idx = 1 : 10
    temp = load(['/tmpscratch/kboerg/aggloGridSearch6/6_01_00046/' sprintf('%0.2u',idx)]);
    bordersGT = [bordersGT; cell2mat(temp.edgesToStore)];
end
edgesGT=find(ismember(graph.borderIdx,bordersGT));
startAgglo = load('/gaba/scratch/mberning/aggloGridSearch/search05_00564.mat');
edgesGTall = [graph.edges(edgesGT,:);startAgglo.axonEdges];
% edgesGTall is bigger than the edges used for the final agglos, because some excluded edges from agglomerateMerge.m made it through (see code there)

result1 = load('/tmpscratch/kboerg/axonQueryAnalysisResult1.mat');
result2 = load('/tmpscratch/kboerg/axonQueryAnalysisResult2.mat');
result1ends = cell2mat(result1.endAgglo(cellfun('length',result1.endAgglo)>0)');
result1ends(result1ends>length(axons)) = [];
a = cellfun(@(x)max([-1,find(ismember(result1ends,x))]),result2.startAgglo);
result2ends = cell2mat(result2.endAgglo(cellfun('length',result2.endAgglo)>0)');
result2ends(result2ends>length(axons)) = [];

edges1 = cellfun(@(x,y)combnk([-1, x(x<=length(axons)) y(y<=length(axons))], 2), result1.startAgglo(result1.idxGood), result1.endAgglo(result1.idxGood), 'uni', 0);
edges2 = cellfun(@(x,y)combnk([-1, x(x<=length(axons)) y(y<=length(axons))], 2), result2.startAgglo(result2.idxGood), result2.endAgglo(result2.idxGood), 'uni', 0);
edges=[edges1; edges2];
edges = cell2mat(edges);
edges(edges(:,1) == edges(:,2),:) = [];
edges(any(edges==-1,2),:)=[];
eqClassCC = Graph.findConnectedComponents(edges, true, true);

eqClassCCfull = [eqClassCC; num2cell(setdiff(1 : length(axons), cell2mat(eqClassCC)))'];
biggestEC =eqClassCC{1000};
biggestECids = cell2mat(axons(biggestEC));
nodes = [];
edges = [];
lookup1 = [];
lookup2 = [];
for idx = 1 : length(biggestEC)
    idx
    lookup=sparse(ones(1,length(axons{biggestEC(idx)})), axons{biggestEC(idx)}, 1:length(axons{biggestEC(idx)}));
    
    edges=[edges;full(lookup(edgesGTall(all(ismember(edgesGTall,axons{biggestEC(idx)}),2),:)))];
    nodes=[nodes;segmentMeta.point(:,axons{biggestEC(idx)})'];
    lookup1 = [lookup1; repmat(biggestEC(idx), length(axons{biggestEC(idx)}), 1)];
    lookup2 = [lookup2; axons{biggestEC(idx)}];
end
edges2=edges;
nodes2=nodes;
for idx = 1 : length(result1.startAgglo)
    idx
    if ~isempty(result1.startAgglo{idx}) && result1.startAgglo{idx}<=length(axons) && ismember(result1.startAgglo{idx},biggestEC)
        tempids =[result1.ff.segIds{idx},result1.ff.neighbours{idx}];
        hits=find(any(ismember(tempids, cell2mat(axons(biggestEC))),2));
        for idx2= 1:length(hits)
            nodes2=[nodes2;result1.ff.nodes{idx}(hits(idx2),:)];
            if idx2>1
                edges2=[edges2;size(nodes2,1),size(nodes2,1)-1];
            end
            edges2=[edges2;size(nodes2,1),find(ismember(lookup2,tempids(hits(idx2),:)),1)];
            
        end
    end
end
