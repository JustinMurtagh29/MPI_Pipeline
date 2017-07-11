% ... for cutting down graph
borderSizeDendrites = 300; % border below this size are excluded from dendrite subgraph during inital agglomeration
segmentSizeDendrites = 500; % segment below ...
borderSizeAxons = 50; % border below this size are excluded from axon subgraph during inital agglomeration
segmentSizeAxons = 100; % segment below ...

% Restrictions of graph based on segment class proability
axonProbThreshold = 0.5;;
dendriteProbThreshold = 0.3;
% Threshold on neurite continuity probability for CC (and final agglomerate size in voxels)
probThresholdDendrite = 0.98;
sizeThresholdDendrite = 1e5;
probThresholdAxon = 0.93;
sizeThresholdAxon = 1e5;

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
connectEM.agglomeration( ... 
        borderSizeDendrites, segmentSizeDendrites, borderSizeAxons, segmentSizeAxons, ...
        axonProbThreshold, dendriteProbThreshold, spineProbThreshold, ...
        probThresholdDendrite, sizeThresholdDendrite, probThresholdAxon, sizeThresholdAxon, ...
        erProbThreshold, ...
        dendriteProbSpines, probThresholdSpines, maxStepsSpines, ...
        '/gaba/scratch/mberning/aggloSearch/00024.mat', graph);
