% generate segmentation from SVM predictions on the whole dataset
rootDir = '/tmpscratch/sahilloo/data/H2_3_v2_U1_SubI/pipelineRun_mr2e_wsmrnet/';
load(fullfile(rootDir, 'allParameter.mat'));

% init wkw dataset, if needed
segSVM.root = '/tmpscratch/sahilloo/data/H2_3_v2_U1_SubI/pipelineRun_mr2e_wsmrnet/SVMpredictions/1/';
segSVM.backend = 'wkwrap';
if ~exist(segSVM.root,'dir')
    mkdir(segSVM.root);
end

wkwInit('new', segSVM.root, 32, 32, 'uint32', 1);

% collect parameters
inputCell = cell(numel(p.local), 1);
for curIdx = 1:numel(p.local)
   [curI, curJ, curK] = ind2sub(size(p.local), curIdx);
    inputCell{curIdx} = {curI, curJ, curK};
end

job = Cluster.startJob( ...
        @H2_3_v2_SubI.SVM.getSegmentationBox, inputCell, ...
        'sharedInputs', {p,segSVM}, ...
        'name', mfilename(), ...
        'taskGroupSize', 10, ...
        'cluster', {'memory', 12});

Cluster.waitForJob(job);

