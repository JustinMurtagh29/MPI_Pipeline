p_col =zeros(0,3);
p_idx =[];
for idx =  1:3000

    if exist([p.local(idx).saveFolder 'segmentAggloPredictions.mat'])
        idx
        pthis = load([p.local(idx).saveFolder 'segmentAggloPredictions.mat']);
        p_idx = [p_idx; pthis.segId];
        p_col = [p_col; pthis.probs];
    end
end
oldPredictions = load('/gaba/u/mberning/results/pipeline/20170217_ROI/segmentAggloPredictionsOld.mat');
newPredictions.segId = unique(cat(1, p_idx, oldPredictions.segId));
[idxA, idxB] = ismember(newPredictions.segId, oldPredictions.segId);
newPredictions.probs(find(idxA),:) = oldPredictions.probs(idxB(idxB>0),:);
[idxA, idxB] = ismember(newPredictions.segId, p_idx);
newPredictions.probs(find(idxA),:) = p_col(idxB(idxB>0),:);
newPredictions.probsMulti = newPredictions.probs ./  repmat(sum(newPredictions.probs,2), 1,3);
Util.saveStruct('/gaba/u/mberning/results/pipeline/20170217_ROI/segmentAggloPredictionsHybrid.mat', newPredictions);

for idx = 1:9
    if  ~isempty(foundlings2{idx})
        foundlings2{idx}=unique(foundlings2{idx});
        disp(idx)
        [(1: sum(ismember(p_idx, foundlings2{idx})))', p_col(ismember(p_idx, foundlings2{idx}), :)]
    end
end
segmentMeta.axonProb(nonull(ygrid5{564}.axon1.ids{1}))

% ... for cutting down graph
borderSizeDendrites = 300; % border below this size are excluded from dendrite subgraph during inital agglomeration
segmentSizeDendrites = 500; % segment below ...
borderSizeAxons = 0; % border below this size are excluded from axon subgraph during inital agglomeration
segmentSizeAxons = 0; % segment below ...

% Restrictions of graph based on segment class proability
axonProbThreshold = 0;
dendriteProbThreshold = 0.3;
% Threshold on neurite continuity probability for CC (and final agglomerate size in voxels)
probThresholdDendrite = 0.98;
sizeThresholdDendrite = 1e5;
probThresholdAxon = 0.96;
sizeThresholdAxon = 0;

% ... for ER reassignment
erProbThreshold = 5;

% NOTE: Parameters currently have no effect, commented routines for first tests
% For deciding in which segments to seed spine attachment greedy search
spineProbThreshold = 0.5;
% ... for spine attachment
dendriteProbSpines = 0.05;
probThresholdSpines = 0.05;
maxStepsSpines = 10;

% Perform agglomeration once first to check whether it works the same way
connectEM.agglomeration(         borderSizeDendrites, segmentSizeDendrites, borderSizeAxons, segmentSizeAxons,     axonProbThreshold, dendriteProbThreshold, spineProbThreshold,         probThresholdDendrite, sizeThresholdDendrite, probThresholdAxon, sizeThresholdAxon,     erProbThreshold,         dendriteProbSpines, probThresholdSpines, maxStepsSpines,    '/gaba/scratch/kboerg/aggloSearch/axonFN96.mat', graph, struct('calculateMetrics', false,'skipDendrites', true));
axonFN99 = load('/gaba/scratch/kboerg/aggloSearch/axonFN99.mat');
axonFN98 = load('/gaba/scratch/kboerg/aggloSearch/axonFN98.mat');
axonFN97 = load('/gaba/scratch/kboerg/aggloSearch/axonFN97.mat');
for idx = 98
    axonFN{idx} = load(['/gaba/scratch/kboerg/aggloSearch/axonFN' num2str(idx) '.mat']);
    foundAgglos{idx} = find(cellfun(@(x)any(ismember(x, allids)), axonFN{idx}.axonsFinal));
    pvalues{idx} = cellfun(@(x)max(segmentMeta.axonProb(x)), axonFN{idx}.axonsFinal(foundAgglos{idx}));

end
for idx = 98
    useAggloCounter{idx} = zeros(1, length(foundAgglos{idx}));
    for idx2 = 1 : length(allids)
        idx2
        usedAgglo = cellfun(@(x)any(ismember(x, allids(idx2))), axonFN{idx}.axonsFinal(foundAgglos{idx}));
        if segmentMeta.axonProb(allids(idx2)) < 0.5
            useAggloCounter{idx}(usedAgglo)=useAggloCounter{idx}(usedAgglo)+1;
        end
        p_col{idx}(idx2) = max([segmentMeta.axonProb(allids(idx2)), pvalues{idx}(usedAgglo)]);
    end
end
sum(cellfun('length', axonFN{98}.axonsFinal)<300)
toswitch = false(size(axonFN{98}.axonsFinal));
segmentMeta2 = segmentMeta;
for idx = 1 : length(axonFN{98}.axonsFinal)
    if ~mod(idx, 1000)
        idx
    end
    toswitch(idx) = length(axonFN{98}.axonsFinal{idx}) < 300 && max(segmentMeta.axonProb(axonFN{98}.axonsFinal{idx})) > 0.9;
    if toswitch(idx)
        segmentMeta2.axonProb(axonFN{98}.axonsFinal{idx}) = 17;
    end
end


connectEM.agglomerationModify(gridAgglo_05{end}, '/gaba/scratch/kboerg/aggloSearch/axonFNrescue.mat', graph, struct('segmentMeta', segmentMeta2, 'skipDendrites', true))
