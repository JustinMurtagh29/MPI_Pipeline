% Old parameters for reference

% ... for cutting down graph
borderSizeDendrites = 100; % border below this size are excluded from dendrite subgraph during inital agglomeration
segmentSizeDendrites = 1000; % segment below ...
borderSizeAxons = 100; % border below this size are excluded from axon subgraph during inital agglomeration
segmentSizeAxons = 1000; % segment below ...

% Restrictions of graph based on segment class proability
axonProbThreshold = 0.9;  
dendriteProbThreshold = 0.5;
% For deciding in which segments to seed spine attachment greedy search
spineProbThreshold = 0.5;

% Threshold on neurite continuity probability (and final agglomerate size in voxels) 
probThresholdDendrite = 0.99;
sizeThresholdDendrite = 1e6;
probThresholdAxon = 0.90;
sizeThresholdAxon = 1e6;

% ... for ER reassignment
erProbThreshold = 2;

% ... for spine attachment
dendriteProbSpines = 0.05;
probThresholdSpines = 0.05;
maxStepsSpines = 10;

% Profile
profile on;

% Perform agglomeration once first to check whether it works the same way
[dendritesFinalWithSpines, dendritesFinal, axonsFinal] = connectEM.agglomeration( ... 
        borderSizeDendrites, segmentSizeDendrites, borderSizeAxons, segmentSizeAxons, ...
        axonProbThreshold, dendriteProbThreshold, spineProbThreshold, ...
        probThresholdDendrite, sizeThresholdDendrite, probThresholdAxon, sizeThresholdAxon, ...
        erProbThreshold, ...
        dendriteProbSpines, probThresholdSpines, maxStepsSpines, ...
        '/gaba/scratch/mberning/profile/test.mat');


profile off;
profsave(profile('info'), '/gaba/scratch/mberning/profile/');

