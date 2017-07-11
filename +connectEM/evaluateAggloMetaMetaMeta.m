yrow = connectEM.evaluateAggloMetaMeta(graph, agglorowall.axonsFinal, agglorowall.dendritesFinal, 'run000all', segmentMeta);
sum(segmentMeta.voxelCount(nonull(axon_13.ids))<100)
nonull=@(x)x(x~=0);
newPredictions.probsMulti(ismember(newPredictions.segId,agglo10.dendritesFinal{3572}),2)
y20.axon.recall_col{2}
numel(nonull(unique(axon_13.ids)))
heuristics.mapping = {heuristics.vesselIdx heuristics.endoIdx heuristics.nucleiIdx ...
    heuristics.addedNucleiIdx heuristics.myelinIdx heuristics.addedMyelinIdx};
heuristics.mapping = cellfun(@find, heuristics.mapping, 'uni', 0);
sum(heuristics.myelinIdx(unique(nonull(axon_13.ids))))
sum(heuristics.vesselIdx(unique(nonull(axon_13.ids))))
sum(heuristics.endoIdx(unique(nonull(axon_13.ids))))
sum(heuristics.nucleiIdx(unique(nonull(axon_13.ids))))
sum(heuristics.addedNucleiIdx(unique(nonull(axon_13.ids))))
sum(heuristics.addedMyelinIdx(unique(nonull(axon_13.ids))))
cellfun(@(x)max(x), graph.neighProb(nonull(unique(axon_13.ids))))
setdiff(nonull(axon_13.ids), cell2mat(agglo20.axonsFinal(y20.axon.foundAgglomerates_col{2})))
cellfun(@(x)max(x), graph.neighProb(ans)



mia = setdiff(setdiff(nonull(axon_13.ids), cell2mat(agglo20.axonsFinal(y20.axon.foundAgglomerates_col{2})))));
find(cellfun(@(x)any(ismember(x, mia)),agglo20.axonsFinal))
newPredictions.probsMulti(ismember(newPredictions.segId,nonull(unique(mh_axon.ids))),:)
find(ismember(newPredictions.segId,nonull(unique(mh_axon.ids))))
ismember(agglo10.dendritesFinal
newPredictions.probsMulti(ismember(newPredictions.segId,nonull(unique(mh_axon.ids)))),:)
newPredictions.segId(ismember(newPredictions.segId,nonull(unique(mh_axon.ids)))))
find(cellfun(@(x)any(ismember(10777414, x)), agglo10.dendritesFinal))
sum(newPredictions.probsMulti(ismember(newPredictions.segId,agglo10.dendritesFinal{3572}),2)>0.9)
for idx = 751:900
    idx
    agglo{idx} = load(['/gaba/scratch/mberning/aggloGridSearch/search03_' sprintf('%0.5u', idx) '.mat'], 'metrics','borderSizeDendrites','segmentSizeDendrites',    'borderSizeAxons',          'segmentSizeAxons',         'axonProbThreshold',      'dendriteProbThreshold',       'spineProbThreshold',     'probThresholdDendrite',     'sizeThresholdDendrite',         'probThresholdAxon',         'sizeThresholdAxon',           'erProbThreshold',        'dendriteProbSpines',       'probThresholdSpines',            'maxStepsSpines');
    %y{idx} = connectEM.evaluateAggloMetaMeta(graph, agglo{idx}.axonsFinal, agglo{idx}.dendritesFinal, ['rungridnew' num2str(idx)], segmentMeta); ;
end
for idx = 751:900
    idx
    if exist(['/gaba/scratch/mberning/aggloGridSearch/search03_' sprintf('%0.5u', idx) '.mat'])
        agglo{idx} = load(['/gaba/scratch/mberning/aggloGridSearch/search03_' sprintf('%0.5u', idx) '.mat'], 'metrics','borderSizeDendrites','segmentSizeDendrites',    'borderSizeAxons',          'segmentSizeAxons',         'axonProbThreshold',      'dendriteProbThreshold',       'spineProbThreshold',     'probThresholdDendrite',     'sizeThresholdDendrite',         'probThresholdAxon',         'sizeThresholdAxon',           'erProbThreshold',        'dendriteProbSpines',       'probThresholdSpines',            'maxStepsSpines');
    else

        agglo{idx} = load(['/gaba/scratch/mberning/aggloGridSearch/search02_' sprintf('%0.5u', idx) '.mat'], 'metrics','borderSizeDendrites','segmentSizeDendrites',    'borderSizeAxons',          'segmentSizeAxons',         'axonProbThreshold',      'dendriteProbThreshold',       'spineProbThreshold',     'probThresholdDendrite',     'sizeThresholdDendrite',         'probThresholdAxon',         'sizeThresholdAxon',           'erProbThreshold',        'dendriteProbSpines',       'probThresholdSpines',            'maxStepsSpines', 'dendritesFinal');
        y{idx} = connectEM.evaluateAggloMetaMeta(graph, agglo{idx}.dendritesFinal, agglo{idx}.dendritesFinal, ['rungridnew' num2str(idx)], segmentMeta); ;
        agglo{idx}.metrics = y{idx};
    end
end

find(cellfun(@(x)x.axonProbThreshold==0.3&&x.probThresholdAxon==0.92&&borderSizeAxons==0))
rmpath ('/gaba/u/kboerg/code/RESCOPaaS/auxiliary/nml')
    rmpath('/gaba/u/kboerg/code/RESCOPaaS/auxiliary')
y676 = connectEM.evaluateAggloMetaMeta(graph, z676.axonsFinal, z676.dendritesFinal, 'run00676', segmentMeta);


find(ismember(segmentMeta.point', [3282, 4201, 1808], 'rows'))
feval(@(x)[graph.edges(x,:), graph.prob(x)*1000], find(any(ismember(graph.edges,14655198),2)))
z676.axonsFinal{y676.axon.foundAgglomerates_col{2}(33)}
temp.graph = graph;
temp.graph.edges = temp.graph.edges(temp.graph.prob>=0.9, :);
yrow = connectEM.evaluateAggloMetaMeta(temp.graph, agglorowall.axonsFinal, agglorowall.dendritesFinal, 'run000all4', segmentMeta); 
agglo21 = load('/gaba/scratch/kboerg/aggloSearch/00021.mat');
y21 = connectEM.evaluateAggloMetaMeta(temp.graph, agglo21.axonsFinal, agglo21.dendritesFinal, 'run00021', segmentMeta);
temp.graph2.edges = agglo22.axonEdges;
temp.graph2.neighbours = graph.neighbours;
y578 = connectEM.evaluateAggloMetaMeta(graph, agglogrid578.axonsFinal, agglogrid578.dendritesFinal, 'run000grid578', segmentMeta); 
y22 = connectEM.evaluateAggloMetaMeta(temp.graph2, agglo22.axonsFinal, agglo22.dendritesFinal, 'run00022b', segmentMeta); 
sum(segmentMeta.voxelCount(agglogrid578.axonsFinal{y578.axon.foundAgglomerates_col{1}(21)}))
agglosplit = load('/gaba/scratch/kboerg/aggloSearch/split.mat');
p

agglogrid2{878} = load('/gaba/scratch/mberning/aggloGridSearch/search02_00878.mat');
globalSynScores = load('/gaba/u/mberning/results/pipeline/20170217_ROI/globalSynScores');
segmentMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta');
graph2 = load('/gaba/u/mberning/results/pipeline/20170217_ROI/graph.mat', 'prob');

ygrid2{878} = connectEM.evaluateAggloMetaMeta(graph, agglogrid2{878}.axonsFinal, agglogrid2{878}.dendritesFinal, 'rungrid2_878', segmentMeta); 
tempbouton = graph.edges(globalSynScores.edgeIdx, :);
tempbouton2 = tempbouton(globalSynScores.synScores>0);
temprandomagglos = agglogrid2{878}.axonsFinal(randperm(length(agglogrid2{878}.axonsFinal), 100));
tempa = cellfun(@(x)length(intersect(x, tempbouton2)), temprandomagglos, 'uni', 0);
segmentMeta2 = segmentMeta;
segmentMeta2.point = segmentMeta2.point';
mkdir('/gaba/scratch/kboerg/throwback_tuesday/');
connectEM.skeletonFromAgglo(graph.edges, segmentMeta2, temprandomagglos(cell2mat(tempa)>0), 'throwback_tuesday', '/gaba/scratch/kboerg/throwback_tuesday/')
tempidx =[];
for idx = 1 : 206770
temp = segmentMeta.point(:, gridAgglo_03{514}.dendritesFinal{idx})';
temp2 = max(temp)-min(temp);
if temp2(1)>4000 && temp2(2)<1000 && temp2(3)<500
    tempidx(end+1) = idx;
end
end
connectEM.skeletonFromAgglo(graph2.edges, segmentMeta2, gridAgglo_03{514}.dendritesFinal(tempidx) , 'apicals', '/gaba/scratch/kboerg/apicals/')
gridAgglo_05{364} = load('/gaba/scratch/mberning/aggloGridSearch/search05_00364.mat');
ygrid5{364} = connectEM.evaluateAggloMetaMeta(graph, gridAgglo_05{364}.axonsFinal, gridAgglo_05{364}.dendritesFinal, 'rungrid5_364', segmentMeta); 

tempa = cellfun(@(x)sum(segmentMeta.axonProb(nonull(x))<0.5),ygrid5{564}.axon1.ids);
tempb = cellfun(@(x)sum(excludedSegmentIdx(nonull(x))),ygrid5{564}.axon1.ids);
tempc = cellfun(@(x)cell2mat(gridAgglo_05{564}.axonsFinal(x)), ygrid5{564}.axon1.foundAgglomerates_col, 'uni', 0);
tempd = cellfun(@(x,y)length(intersect(y, nonull(x))),ygrid5{564}.axon1.ids,tempc);
tempe = cellfun(@(x)length(nonull(x)),ygrid5{564}.axon1.ids);
for idx = 1 : 10
    disp(sprintf('axon %u: %u segments, %u picked up, %u axonexcluded, %u from bad cubes, %u singular and excluded for other reasons',idx,tempe(idx), tempd(idx), tempa(idx), tempb(idx), tempe(idx)-tempd(idx)-tempa(idx)-tempb(idx)))
end
sum(ygrid5{564}.axon1.ids{1}==0)
tempaa = cellfun(@(x)segmentMeta.point(:, feval(@(y)y(segmentMeta.axonProb(y)<0.5), nonull(x)))',ygrid5{564}.axon1.ids, 'uni', 0);
tempdd = cellfun(@(x,y)segmentMeta.point(:, setdiff(nonull(x), y))',ygrid5{564}.axon1.ids,tempc, 'uni', 0);
aggloposthoc=load('/gaba/scratch/kboerg/aggloSearch/aggloPostHoc');
yposthoc = connectEM.evaluateAggloMetaMeta(graph, aggloposthoc.axonsFinal, aggloposthoc.dendritesFinal, 'aggloposthoc', segmentMeta); 
y_zz = connectEM.evaluateAggloMetaMeta(graph, zz.axonsFinal, zz.dendritesFinal, 'cycle3', segmentMeta);
liste=dir(['/gaba/scratch/kboerg/eval_agglo/run564again2/axons1/3_*.nml'])
mkdir('/gaba/scratch/kboerg/amira2/')
for idx = 1 :length(liste)
    idx
    galleryCortex(p, '/gaba/scratch/kboerg/eval_agglo/run564again2/axons1/', liste(idx).name, '/gaba/scratch/kboerg/amira3/');
end
cycleAgglo_01{3} = load(['/gaba/scratch/kboerg/directCycle/cycle' num2str(3, '%0.3u')]);
cycleAgglo_01{4} = load(['/gaba/scratch/kboerg/directCycle/cycle' num2str(4, '%0.3u')]);
y_cycleAgglo_01{6} = connectEM.evaluateAggloMetaMeta(graph, cycleAgglo_01{6}.axonsFinal, cycleAgglo_01{6}.dendritesFinal, 'cycleAgglo_01_06', segmentMeta);
y_cycleAgglo_01{7} = connectEM.evaluateAggloMetaMeta(graph, cycleAgglo_01{7}.axonsFinal, cycleAgglo_01{7}.dendritesFinal, 'cycleAgglo_01_07', segmentMeta);


directions_cycle{5} = load('/gaba/scratch/kboerg/directCycle/temp005/directions');

% probability
% check border prob
feval(@(x,y)graph.neighProb{x}(graph.neighbours{x} == y),  2787403,   2787768)
% check border size
feval(@(x,y)borderMeta.borderSize(graph.neighBorderIdx{x}(graph.neighbours{x} == y)), 2787403,   2787768)
% find edge id
feval(@(x,y)graph.neighBorderIdx{x}(graph.neighbours{x} == y), 11134065,   11134338)

% check latent score
directions.latent(2787617)
% check member of agglo
find(cellfun(@(x)any(ismember(x, 11134065)), axons))
agglos{530826}
%check agglo size
directions.agglomerationSize(2787617)
%check border score
feval(@(x)directions.scores(x), ismember(directions.edges, [11134338,11134065], 'rows'))
% in directions idx
find(ismember(directions.edges, [11134338,11134065], 'rows'))
find(ismember(graph.edges, [9945754,11134065], 'rows'))



feval(@(x)directions_cycle{5}.directions.scores(x), ismember(directions_cycle{5}.directions.edges, fliplr([603917, 604321]), 'rows'))
forceKeepEdges{3} = load('/gaba/scratch/kboerg/directCycle/forceKeepEdges_003');
find(ismember(graph.edges(forceKeepEdges{3}.forceKeepEdgesStore,:),[11134065, 11134338],'rows'))
segmentMeta.voxelCount(4162169)
segmentMeta.axonProb(   13842294)
forcingNum(idx) = connectEM.agglomerationPostHocTwo(options, [topfolder, 'cycle', num2str(idx + 1, '%0.3u')], graph, borderMeta, segmentMeta, directions_cycle{5}.directions, cycleAgglo_01{6},idx,topfolder);
forceCorrespondences(71337821)
find(ismember(directions.edges, [604321, 603917], 'rows'))
find(ismember(optional.forceKeepEdges, [9945754,11134065], 'rows'))
find(optional.forceKeepEdges== 71337821)
find(ismember(segmentMeta.point, [2406, 4932, 669]+1, 'rows'))
yAM = connectEM.evaluateAggloMetaMeta(graph, aggloAM.axonsMerged, [], 'aggloAM', segmentMeta); 
