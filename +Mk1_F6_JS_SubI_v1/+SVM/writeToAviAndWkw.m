% write a test bbox prediction to webknossos
% These predictions were generated from mr_1g_cont5 cnn with python
% This script generates avi movie and webknossos wkw segmentation

addpath(genpath('/u/sahilloo/repos/Benedikt/'));
datasetName = 'Mk1_F6_JS_SubI_v1';
svmPredNetwork = 'mr_2e_cont4';
boxId = 1;

rootDir = fullfile('/tmpscratch/sahilloo/data/', datasetName, 'pipelineRun_mr2e_wsmrnet/');
m = load(fullfile(rootDir, 'allParameter.mat'));
p = m.p;

outputDirWkw = fullfile('/tmpscratch/sahilloo/data/', datasetName, ['svm_' svmPredNetwork '_wkw']);
mkdir(outputDirWkw)

outputDir = fullfile('/tmpscratch/sahilloo/data/', datasetName, ['svm_' svmPredNetwork '_mat']);
mkdir(outputDir)

Util.log('load python svm predictions')
m = load(fullfile('/tmpscratch/sahilloo/data/', datasetName, ['svm_' svmPredNetwork '_slurm'], [num2str(boxId) '.mat']));
pred = m.pred;
bbox_wk = [m.offset, m.shape];
bbox = Util.convertWebknossosToMatlabBbox(bbox_wk);

predOut = zeros([bbox_wk(4:end), 3]);
predOut(:,:,:,2) = pred(:,:,:,3);
predOut(:,:,:,3) = pred(:,:,:,4);

% obtain raw
raw = Seg.IO.loadRaw(p, bbox);

Util.log('obtain seg')
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

Util.log('make avi file')
Visualization.movieMakerSeg(raw,seg,fullfile(outputDir,['seg_svm_' num2str(boxId,'%05d') '.avi']))

% write to WKW
segSVM.root = fullfile(outputDirWkw,'1');
segSVM.backend = 'wkwrap';
if ~exist(segSVM.root,'dir')
    mkdir(segSVM.root);
end
wkwInit('new',segSVM.root,32, 32, 'uint32', 1);
saveSegDataGlobal(segSVM, double(bbox(:,1)'), seg); 

