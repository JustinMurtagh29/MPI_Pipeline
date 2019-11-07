% write a test bbox prediction to webknossos
% These predictions were generated from  cnn17 with python

addpath(genpath('/u/sahilloo/repos/Benedikt/'));

rootDir = '/tmpscratch/sahilloo/data/H2_3_v2_U1_SubI/pipelineRun_mr2e_wsmrnet/';
m = load(fullfile(rootDir, 'allParameter.mat'));
p = m.p;

% write to WKW
segSVM.root = fullfile(rootDir,'SVMPredictionsWkwTestCNN17/1/');
segSVM.backend = 'wkwrap';
if ~exist(segSVM.root,'dir')
    mkdir(segSVM.root);
end

wkwInit('new',segSVM.root,32, 32, 'uint32', 1);

Util.log('load matlab svm predictions')
idxDebug = [22    44    66    88   110   132   154   176   198   220];
for index=1:numel(idxDebug)
    m = load(fullfile(p.local(idxDebug(index)).saveFolder, 'svmPredictions.mat'),...
                 'pred17', 'pred17_2', 'pred17_4', 'bbox');
    bbox = m.bbox;
    m1_pred = m.pred17_2; % use for mitos
    m2_pred = m.pred17_4; % use for vesicles
    
    bbox_wk = Util.convertMatlabToWKBbox(bbox)
    
    predOut = zeros([bbox_wk(4:end), 3]);
    predOut(:,:,:,2) = m2_pred(:,:,:,3);
    predOut(:,:,:,3) = m1_pred(:,:,:,4);
    
    % obtain raw
    raw = Seg.IO.loadRaw(p, bbox);
    
    %Util.log('obtain seg')
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
    
    %Util.log('make avi file')
    Visualization.movieMakerSeg(raw,seg,fullfile(rootDir,['seg_svm_cnn17_' num2str(index) '.avi']))
    
    saveSegDataGlobal(segSVM, double(bbox(:,1)'), seg); 
end
