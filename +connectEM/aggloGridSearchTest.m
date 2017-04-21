% ... for cutting down graph
borderSizeDendrites = 100; % border below this size are excluded from dendrite subgraph during inital agglomeration
segmentSizeDendrites = 100; % segment below ...
borderSizeAxons = 100; % border below this size are excluded from axon subgraph during inital agglomeration
segmentSizeAxons = 100; % segment below ...

% Restrictions of graph based on segment class proability
axonProbThreshold = 0.3;
dendriteProbThreshold = 0;
% For deciding in which segments to seed spine attachment greedy search
spineProbThreshold = 0.5;

% Threshold on neurite continuity probability for CC (and final agglomerate size in voxels) 
probThresholdDendrite = 0.99;
sizeThresholdDendrite = 100;
probThresholdAxon = 0.90;
sizeThresholdAxon = 100;

% NOTE: Parameters currently have no effect, commented routines for first tests
% ... for ER reassignment
erProbThreshold = 2;
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
        '/gaba/scratch/mberning/aggloSearch/00020.mat');

