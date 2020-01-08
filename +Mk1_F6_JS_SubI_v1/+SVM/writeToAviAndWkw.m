% write a test bbox prediction to webknossos
% These predictions were generated from mr_1g_cont5 cnn with python
% This script generates avi movie and webknossos wkw segmentation

addpath(genpath('/u/sahilloo/repos/Benedikt/'));
datasetName = 'Mk1_F6_JS_SubI_v1';
svmPredNetwork = 'mr_2e_cont4';
dateStamp =  datestr(clock,30);

aviFlag = false;
rootDir = fullfile('/tmpscratch/sahilloo/data/', datasetName, 'pipelineRun_mr2e_wsmrnet/');
m = load(fullfile(rootDir, 'allParameter.mat'));
p = m.p;

outputDirWkw = fullfile('/tmpscratch/sahilloo/data/', datasetName, ['svm_' svmPredNetwork '_wkw'], dateStamp);
mkdir(outputDirWkw)

outputDir = fullfile('/tmpscratch/sahilloo/data/', datasetName, ['svm_' svmPredNetwork '_mat'], dateStamp);
mkdir(outputDir)

% write to WKW
segSVM.root = fullfile(outputDirWkw,'1');
segSVM.backend = 'wkwrap';
if ~exist(segSVM.root,'dir')
    mkdir(segSVM.root);
end
wkwInit('new',segSVM.root,32, 32, 'uint32', 1);

boxIds = 1:prod(p.tiles);
tic;
endStep = numel(boxIds);
for boxId = 1:endStep
    m = load(fullfile('/tmpscratch/sahilloo/data/', datasetName, ['svm_' svmPredNetwork '_slurm'], [num2str(boxId) '.mat']));
    pred = m.pred;
    bbox_wk = [m.offset, m.shape];
    bbox = Util.convertWebknossosToMatlabBbox(bbox_wk);
    
    predOut = zeros([bbox_wk(4:end), 3]); % segId 1 -> PSD, here it is empty
    predOut(:,:,:,2) = pred(:,:,:,3); % segId 2 -> VC
    predOut(:,:,:,3) = pred(:,:,:,4); % segId 3 -> Mito
    
    Util.log('obtain seg')
    t = 0.5;
    minArea = 3000;
    nhood = [5,5,3];
    [seg, synCom, vcCom, miCom] = ...
        Paper.SynEM.MethComp.preprocessSVMForAnnotation(predOut, t, minArea, nhood);
    
    % further delete small components
    for i = 1:3
        tmp = seg == i;
        tmp = regionprops(tmp > 0, 'Area', 'PixelIdxList');
        tmp([tmp.Area] >= minArea) = [];
        idx = cell2mat({tmp.PixelIdxList}');
        seg(idx) = 0;
    end
    if aviFlag
        Util.log('make avi file')
        % obtain raw
        raw = Seg.IO.loadRaw(p, bbox);
        Visualization.movieMakerSeg(raw,seg,fullfile(outputDir,['seg_svm_' num2str(boxId,'%05d') '.avi']))
    end
    saveSegDataGlobal(segSVM, double(bbox(:,1)'), seg);
    Util.save([p.local(boxId).saveFolder 'svmSegData.mat'], seg, synCom, vcCom, miCom);
    Util.progressBar(boxId, endStep);
end
% create all resolutions for segmentation
%thisBBox = [1, 1, 1; (ceil(p.bbox(:, 2) ./ 1024) .* 1024)']';
%createResolutionPyramid(segSVM, thisBBox, [], true);
