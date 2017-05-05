topfolder = '/gaba/scratch/kboerg/directCycle/';
mkdir(topfolder);

% initialize
forceKeepEdgesStore = [];
save([topfolder, 'forceKeepEdges_000'], 'forceKeepEdgesStore');
clear forceKeepEdgesStore;
copyfile('/gaba/scratch/mberning/aggloGridSearch/search05_00564.mat', [topfolder, 'cycle001.mat']);

% load big structures
borderMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/globalBorder.mat', 'borderSize', 'borderCoM');
segmentMeta = load([p.saveFolder 'segmentMeta.mat'], 'voxelCount', 'point', 'maxSegId', 'cubeIdx');
segmentMeta.point = segmentMeta.point';
segmentMeta = connectEM.addSegmentClassInformation(p, segmentMeta);

% set options
options.borderSizeThreshold = 30;
options.sizeThreshold = 100;
options.axonProbThreshold = 0.5000;
options.probThreshold = 0.8000;
options.scoreThreshold = 0.9000;
options.latentThreshold = 0.7000;
options.agglomerationSizeThreshold = 2500;

% iterate
tempfolder = @(k)[topfolder 'temp' , num2str(k, '%0.3u') filesep];
for idx = 1 : 10
    mkdir(tempfolder(idx));
    agglo = load([topfolder, 'cycle', num2str(idx, '%0.3u')]);
    job = connectEM.agglomerateDirectionalitySuperStart([topfolder, 'cycle', num2str(idx, '%0.3u')], 0.5, tempfolder(idx));
    Cluster.waitForJob(job);
    directions = connectEM.agglomerateDirectionalityCollect(tempfolder(idx), segmentMeta, agglo, graph);
    forcingNum(idx) = connectEM.agglomerationPostHocTwo(options, [topfolder, 'cycle', num2str(idx + 1, '%0.3u')], graph, borderMeta, segmentMeta, directions, agglo,idx,topfolder);
end
