function makeSkelForChiasmataDetection(startingidx)
load('/tmpscratch/kboerg/axonsBorderNew');
segmentMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta.mat', 'point');
load('/gaba/scratch/mberning/edgesGTall.mat');
result3 = load('/tmpscratch/kboerg/axonQueryAnalysisResult3.mat');




edges = cellfun(@(x,y)combnk([-1, x(x<=length(axons)) y(y<=length(axons))], length(x(x<=length(axons)))*2), result3.startAgglo(result3.idxGood), result3.endAgglo(result3.idxGood), 'uni', 0);
edges(cellfun(@isempty,edges))=[];
edges = cell2mat(edges);
edges(edges(:,1) == edges(:,2),:) = [];
edges(any(edges==-1,2),:)=[];
eqClassCC = Graph.findConnectedComponents(edges, true, true);

eqClassCCfull = [eqClassCC; num2cell(setdiff(1 : length(axons), cell2mat(eqClassCC)))'];
% iterate over super agglos
numstr = '19';

for idx_agglo = startingidx : 500 : length(eqClassCCfull);
    currentEC =eqClassCCfull{idx_agglo};
    mkdir(['/tmpscratch/kboerg/visX', numstr, '_' num2str(floor(idx_agglo/100)) '/']);
    [nodes2, edges2] = connectEM.makeSkelForChiasmataDetectionSub(currentEC, axons, edgesGTall, {result3}, segmentMeta, false);
    assert(length(Graph.findConnectedComponents(edges2))<=1);
    connectEM.detectChiasmata([],nodes2,edges2,idx_agglo ~= 1,['/tmpscratch/kboerg/visX',numstr ,'_' num2str(floor(idx_agglo/100)) '/visX', numstr, '_' num2str(idx_agglo) '/'])
end