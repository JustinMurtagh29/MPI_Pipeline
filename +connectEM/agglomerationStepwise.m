function agglomerationStepwise(graph, segmentMeta)
borderSizeDendrites = 500; % border below this size are excluded from dendrite subgraph during inital agglomeration
segmentSizeDendrites = 500; % segment below ...
dendriteProbThreshold = 0.5; % restrictions of graph based on segment class proability
probThresholdDendrite = 0.99; % threshold on neurite continuity probability for CC
sizeThresholdDendrite = 100; % theshold on final agglomerate size in voxels (single segments not collected!)

% axon related parameter in this run
borderSizeAxons = 0; % border below this size are excluded from axon subgraph during inital agglomeration
segmentSizeAxons = 0; % segment below ...
axonProbThreshold = 0.3; % restrictions of graph based on segment class proability
probThresholdAxon = 0.98; % threshold on neurite continuity probability for CC
sizeThresholdAxon = 0; % theshold on final agglomerate size in voxels (single segments not collected!)

% NOTE: Parameters currently have no effect, commented routines for first tests
% ... for ER reassignment
erProbThreshold = 2;
% ... for spine attachment
spineProbThreshold = 0.5;
dendriteProbSpines = 0.05;
probThresholdSpines = 0.05;
maxStepsSpines = 10;
segmentMeta2=segmentMeta;
for idx = 98:-2:86
    probThresholdAxon = idx/100;
    connectEM.agglomeration(borderSizeDendrites, segmentSizeDendrites, borderSizeAxons, segmentSizeAxons,axonProbThreshold, dendriteProbThreshold, spineProbThreshold,probThresholdDendrite, sizeThresholdDendrite, probThresholdAxon, sizeThresholdAxon, erProbThreshold, dendriteProbSpines, probThresholdSpines, maxStepsSpines, ['/gaba/scratch/kboerg/aggloSearch/' num2str(idx) '.mat'], graph, segmentMeta2);
    agglorow{idx} = load(['/gaba/scratch/kboerg/aggloSearch/' num2str(idx) '.mat']);
    toobigrow{idx} = cellfun(@length,agglorow{idx}.axonsFinal)>100;
    segmentMeta2.voxelCount(cell2mat(agglo96.axonsFinal(toobigrow{idx}))) = -1;
end
