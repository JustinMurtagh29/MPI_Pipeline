function makeSkelForChiasmataDetection(startingidx)
temp = load('/gaba/scratch/mberning/axonQueryGeneration/beforeQueryGeneration.mat', 'axonsNew');
axons = temp.axonsNew;
clear temp
segmentMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta.mat', 'point');
load('/gaba/scratch/mberning/edgesGTall.mat');
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
last = @(x)x{end};
fifth = @(x)x{5}
getTask = @(x)fifth(strsplit(last(strsplit(x,'/')),'_'));
% iterate over super agglos
for idx_agglo = startingidx : 500 : length(eqClassCCfull)
    currentEC =eqClassCCfull{idx_agglo};
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
    resultCol = {result1, result2};
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
    mkdir(['/tmpscratch/kboerg/visX11_' num2str(floor(idx_agglo/100)) '/']);
    detectChiasmataKMB2([],nodes2,edges2,true,['/tmpscratch/kboerg/visX11_' num2str(floor(idx_agglo/100)) '/visX11_' num2str(idx_agglo) '/'])
end