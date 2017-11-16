function getMitoScoreAndVolume(p)

% Get mito scores
labelMap(1).root = '/tmpscratch/bstaffle/data/2012-09-28_ex145_07x2_ROI2017/SVM/cnn17_2';
labelMap(1).backend = 'wkwrap';
labelMap(1).channel = 4;
labelMap(1).segId = [0.5 1];
labelMap(1).function = @(x) bwareaopen(x,1000); % function that is applied to the label map after thresholdding


inputCell = cell(numel(p.local),1);
for cube=1:numel(p.local)
    bbox = p.local(cube).bboxSmall;
    inputCell{cube} = {bbox};
end
functionH = @connectEM.getSegmentHeuristicsScore;

% Execute and retrieve results
seg = p.seg;
cluster = Cluster.getCluster('-l h_vmem=12G');
job = Cluster.startJob(functionH, inputCell, ...
    'name', 'mitoLookup', 'sharedInputs', {seg labelMap}, ...
    'sharedInputsLocation', [1 2], 'cluster', cluster, ...
    'numOutputs', 3);
Cluster.waitForJob(job);

scores = fetchOutputs(job);

% Reformatting for save file
segIds = cat(1, scores{:,2});
sc = cat(1, scores{:,1});
mitoScore = cat(1, sc{:,1});
sc = cat(1, scores{:,3});
voxelCounts = cat(1, sc{:,1});

Util.save([p.saveFolder 'mitoResult.mat'], segIds, mitoScore, voxelCounts);