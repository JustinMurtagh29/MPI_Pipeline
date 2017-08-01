function makeSkelForChiasmataDetection(startingidx)
temp = load('/gaba/scratch/mberning/axonQueryGeneration/beforeQueryGeneration.mat', 'axonsNew');
axons = temp.axonsNew;
clear temp
segmentMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta.mat', 'point');
load('/gaba/scratch/mberning/edgesGTall.mat');
result1 = load('/tmpscratch/kboerg/axonQueryAnalysisResult1.mat');
result2 = load('/tmpscratch/kboerg/axonQueryAnalysisResult2.mat');
filterDoubles = false;
if filterDoubles
    result1ends = cell2mat(result1.endAgglo(cellfun('length',result1.endAgglo)>0)');
    result1ends(result1ends>length(axons)) = [];
    a = find(cellfun(@(x)max([-1,find(ismember(result1ends,x))]),result2.startAgglo)~=-1);
    result2.endAgglo(a) = [];
    result2.startAgglo(a) = [];
    result2.idxGood(a) = [];
end

resultCol = {result1, result2};

edges1 = cellfun(@(x,y)combnk([-1, x(x<=length(axons)) y(y<=length(axons))], length(x(x<=length(axons)))*2), result1.startAgglo(result1.idxGood), result1.endAgglo(result1.idxGood), 'uni', 0);
edges2 = cellfun(@(x,y)combnk([-1, x(x<=length(axons)) y(y<=length(axons))], length(x(x<=length(axons)))*2), result2.startAgglo(result2.idxGood), result2.endAgglo(result2.idxGood), 'uni', 0);
edges=[edges1; edges2];
edges(cellfun(@isempty,edges))=[];
edges = cell2mat(edges);
edges(edges(:,1) == edges(:,2),:) = [];
edges(any(edges==-1,2),:)=[];
eqClassCC = Graph.findConnectedComponents(edges, true, true);

eqClassCCfull = [eqClassCC; num2cell(setdiff(1 : length(axons), cell2mat(eqClassCC)))'];
% iterate over super agglos
numstr = '17';
for idx_agglo = startingidx : 500 : length(eqClassCCfull);
    currentEC =eqClassCCfull{idx_agglo};
    mkdir(['/tmpscratch/kboerg/visX', numstr, '_' num2str(floor(idx_agglo/100)) '/']);
    [nodes2, edges2] = connectEM.makeSkelForChiasmataDetectionSub(currentEC, axons, edgesGTall, resultCol, segmentMeta, false);
    if idx_agglo == 1
        save('backup','nodes2','edges2');
    end
    assert(length(Graph.findConnectedComponents(edges2))<=1);
    connectEM.temp2([],nodes2,edges2,true,['/tmpscratch/kboerg/visX',numstr ,'_' num2str(floor(idx_agglo/100)) '/visX', numstr, '_' num2str(idx_agglo) '/'])
end