% write a test bbox prediction to webknossos

addpath(genpath('/u/sahilloo/repos/Benedikt/'));

rootDir = '/tmpscratch/sahilloo/data/H2_3_v2_U1_SubI/pipelineRun_mr2e_wsmrnet/';
m = load(fullfile(rootDir, 'allParameter.mat'));
p = m.p;

% load python svm predictions
m = load('/tmpscratch/sahilloo/data/H2_3_v2_U1_SubI/svm_mr_1g_cont5_slurm/old/0.mat');
pred = m.pred;
bbox_wk = [m.offset, m.shape];
bbox = Util.convertWebknossosToMatlabBbox(bbox_wk);

predOut = zeros([bbox_wk(4:end), 3]);
predOut(:,:,:,2) = pred(:,:,:,3);
predOut(:,:,:,3) = pred(:,:,:,4);

% obtain raw
raw = Seg.IO.loadRaw(p, bbox);

% obtain seg
t = 0.5;
minArea = 3000;
nhood = [5,5,3];
seg = SVM.getSVMSegmentation(predOut, t, minArea, nhood);

% further delete small components
for i = 1:3
    tmp = seg == i;
    tmp = regionprops(tmp > 0, 'Area', 'PixelIdxList');
    tmp([tmp.Area] >= minArea) = [];
    idx = cell2mat({tmp.PixelIdxList}');
    seg(idx) = 0;
end

% make avi file
Visualization.movieMakerSeg(raw,seg,fullfile(rootDir,'seg_svm.avi'))


