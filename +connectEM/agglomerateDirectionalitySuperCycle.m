topfolder = '/gaba/scratch/kboerg/directCycle/';
mkdir(topfolder);
copyfile('/gaba/scratch/mberning/aggloGridSearch/search05_00564.mat', [topfolder, 'cycle001.mat']);
borderMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/globalBorder.mat', 'borderSize', 'borderCoM');
segmentMeta = load([p.saveFolder 'segmentMeta.mat'], 'voxelCount', 'point', 'maxSegId', 'cubeIdx');
segmentMeta.point = segmentMeta.point';
segmentMeta = connectEM.addSegmentClassInformation(p, segmentMeta);
options.borderSizeThreshold = 30;
options.sizeThreshold = 100;
options.axonProbThreshold = 0.5000;
options.probThreshold = 0.8000;
options.scoreThreshold = 0.9000;
options.latentThreshold = 0.7000;
options.agglomerationSizeThreshold = 2500;

for idx = 1 : 10
    tempfolder = [topfolder 'temp' , num2str(idx, '%0.3u')];
    mkdir(tempfolder);
    agglo = load([topfolder, 'cycle', num2str(idx, '%0.3u')]);
    job = agglomerateDirectionalitySuperStart([topfolder, 'cycle', num2str(idx, '%0.3u')], 0.5, tempfolder);
    while isempty(job.FinishTime)
        pause(1);
    end
    directions = agglomerateDirectionalityCollect(tempfolder);
    agglomerationPostHocTwo(options, [topfolder, 'cycle', num2str(idx, '%0.3u')], graph, borderMeta, segmentMeta, directions, agglo);
end
