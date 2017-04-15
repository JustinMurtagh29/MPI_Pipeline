% Load parameters
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameter.mat');

job = correspondenceFinderGlobal(p);
Cluster.waitForJob(job);

corrEdges = connectEM.collectGlobalCorrespondences(p);

cc = Graph.findConnectedComponents(corrEdges, true, true);
segMeta = load([p.saveFolder 'segmentMeta.mat'], 'maxSegId');
script = WK.makeMappingScript(segMeta.maxSegId, cc, false);
fileHandle = fopen('/gaba/scratch/mberning/newCorrespondences.txt', 'w+');
fwrite(fileHandle, script);
fclose(fileHandle);

save([p.saveFolder 'correspondencesNew.mat'], 'corrEdges');

